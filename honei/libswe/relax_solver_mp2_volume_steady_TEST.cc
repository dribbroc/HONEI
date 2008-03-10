/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/libswe/relax_solver.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <string>
#include <honei/libswe/scenario.hh>
#include <sys/time.h>
#include <volume.hh>
#include <honei/libmath/quadrature.hh>
using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;
using namespace precision_modes;

template <typename Tag_>
class RelaxSolverMIXEDPRECINNERVolTest :
    public BaseTest
{
    public:
        RelaxSolverMIXEDPRECINNERVolTest(const std::string & type) :
            BaseTest("relax_solver_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            ulint dwidth = 40;
            ulint dheight = 40;
            ulint timesteps = 5000;

            DenseMatrix<float> height(dheight, dwidth, float(5));
            DenseMatrix<float> bottom(dheight, dwidth, float(1));

            Cuboid<float> a(height, float(5), float(5), float(10), 20, 20);
            a.value();

            DenseMatrix<float> u1(dheight, dwidth, float(0));
            DenseMatrix<float> u2(dheight, dwidth, float(0));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<float> u(entries, float(1));
            DenseVector<float> v(entries, float(1));
            DenseVector<float> w(entries, float(1)); 
            DenseVector<float> bx (entries/3, float(0));
            DenseVector<float> by (entries/3, float(0));
            DenseVector<float> c (3,float(5));
            DenseVector<float> d (3,float(5));
            //SCENARIO setup:
            c[0] = 10;
            c[1] = 6;
            c[2] = 11;
            d[0] = 10;
            d[1] = 5;
            d[2] = 11;

            float deltax = 5;
            float deltay = 5;
            float deltat = 5./24.;

            double eps = 10e-6;
            float manning = float(0);
            Scenario<float, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);
            scenario.height = &height;
            scenario.bottom = &bottom;
            scenario.x_veloc = &u1;
            scenario.y_veloc = &u2;
            scenario.u = &u;
            scenario.v = &v;
            scenario.w = &w;
            scenario.bottom_slopes_x = &bx;
            scenario.bottom_slopes_y = &by;
            scenario.c = &c;
            scenario.d = &d;
            scenario.delta_x = deltax;
            scenario.delta_y = deltay;
            scenario.delta_t = deltat;
            scenario.eps = eps;
            scenario.manning_n = manning;

            DenseVector<double> bx_2(bx.size());
            DenseVector<double> by_2(by.size());
            convert(bx_2, bx);
            convert(by_2, by);
            RelaxSolver<Tag_, float, float, double, float, float, source_types::SIMPLE, boundaries::REFLECT, MIXED> relax_solver
                (scenario, bx_2, by_2);
            relax_solver.do_preprocessing();

            DenseMatrix<double> result(height.rows(), height.columns());
#ifdef SOLVER_POSTPROCESSING_VOLUME
            double volumes[timesteps + 1];
            volumes[0] = 0.;
#endif

            for (ulint i = 1; i <= timesteps; ++i)
            {
                relax_solver.solve();
#ifdef SOLVER_POSTPROCESSING_VOLUME
                convert(result, height);
                volumes[i] = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(result, double(0), double(deltax * dwidth), (double)deltax, (double)deltay);
#endif
                cout << "Timestep "<< i <<" / " << timesteps << " finished." <<endl;
            }

            convert(result, height);
            double ana_vol = 0.5 * a.size()* deltax * deltay + (dwidth * deltax * dheight * deltay * 5.);
            std::cout << "Analytical start: " << ana_vol;
            std::cout << "Analytical target: " << ana_vol - 0.5 * a.size()* deltax * deltay<< std::endl<< std::endl;
            double vol =(double)GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(result, double(0), double(deltax * dwidth), (double)deltax, (double)deltay);
            std::cout << "Vol.: " << vol << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(vol, (ana_vol - 0.5 * a.size()* deltax * deltay), 2.);
#ifdef SOLVER_POSTPROCESSING_VOLUME
            std::string filename;
            std::ofstream file;
            filename = "out_relax_vol_mixed2.dat";
            file.open(filename.c_str());
            for(unsigned long i(1); i <= timesteps; ++i)
            {
                file << stringify(i) + " " + stringify(volumes[i]) + " " + stringify(fabs(volumes[i] - (ana_vol - 0.5 * a.size()* deltax * deltay) ))  +"\n";
            }
            file.close();
#endif

        }
};
#ifdef HONEI_SSE
RelaxSolverMIXEDPRECINNERVolTest<tags::CPU::SSE> sse_relax_solver_mp2_test("sse mixedprec variant 2");
#endif
#ifdef HONEI_CELL
RelaxSolverMIXEDPRECINNERVolTest<tags::Cell> cell_relax_solver_mp2_test("CELL mixedprec variant 2");
#endif
