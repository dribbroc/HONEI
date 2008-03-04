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
#include <iostream>
#include <fstream>

using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;
using namespace precision_modes;

template <typename Tag_, typename DataType_>
class RelaxSolverVolumeSteadyTest :
    public BaseTest
{
    public:
        RelaxSolverVolumeSteadyTest(const std::string & type) :
            BaseTest("relax_solver_volume_steady_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            ulint dwidth = 40;
            ulint dheight = 40;
            ulint timesteps = 5000;

            DenseMatrix<DataType_> height(dheight, dwidth, DataType_(5));
            DenseMatrix<DataType_> bottom(dheight, dwidth, DataType_(1));

            Cuboid<DataType_> a(height, DataType_(5), DataType_(5), DataType_(10), 20, 20);
            a.value();

            DenseMatrix<DataType_> u1(dheight, dwidth, DataType_(0));
            DenseMatrix<DataType_> u2(dheight, dwidth, DataType_(0));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<DataType_> u(entries, DataType_(1));
            DenseVector<DataType_> v(entries, DataType_(1));
            DenseVector<DataType_> w(entries, DataType_(1)); 
            DenseVector<DataType_> bx (entries/3, DataType_(0));
            DenseVector<DataType_> by (entries/3, DataType_(0));
            DenseVector<DataType_> c (3,DataType_(5));
            DenseVector<DataType_> d (3,DataType_(5));
            //SCENARIO setup:
            c[0] = 10;
            c[1] = 6;
            c[2] = 11;
            d[0] = 10;
            d[1] = 5;
            d[2] = 11;

            DataType_ deltax = 5;
            DataType_ deltay = 5;
            DataType_ deltat = 5./24.;

            double eps = 10e-6;
            DataType_ manning = DataType_(0);
            Scenario<DataType_, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);
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

            RelaxSolver<Tag_, DataType_, DataType_, DataType_, DataType_, DataType_, source_types::SIMPLE, boundaries::REFLECT, FIXED> relax_solver
                (scenario);
            relax_solver.do_preprocessing();

#ifdef SOLVER_POSTPROCESSING_VOLUME
            DataType_ volumes[timesteps + 1];
            volumes[0] = 0.;
#endif
            for (ulint i = 1; i <= timesteps; ++i)
            {
                relax_solver.solve();
#ifdef SOLVER_POSTPROCESSING_VOLUME
                volumes[i] = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(height, DataType_(0), DataType_(deltax * dwidth), deltax, deltay);
#endif
                cout << "Timestep "<< i <<" / " << timesteps << " finished. " <<endl;
            }
            DataType_ ana_vol = 0.5 * a.size()* deltax * deltay + (dwidth * deltax * dheight * deltay * 5.);
            std::cout << "Analytical start: " << ana_vol;
            std::cout << "Analytical target: " << ana_vol - 0.5 * a.size()* deltax * deltay<< std::endl<< std::endl;
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(height, DataType_(0), DataType_(deltax * dwidth), deltax, deltay);
            std::cout << "Vol.: " << vol << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(vol, (ana_vol - 0.5 * a.size()* deltax * deltay), 0.3);

#ifdef HONEI_POSTPROCESSING_VOLUME
            std::string filename;
            std::ofstream file;
            filename = "out_relax_vol_fixed.dat";
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
RelaxSolverVolumeSteadyTest<tags::CPU::SSE, double> sse_relax_solver_vs_test_double("VS sse double");
#endif

#ifdef HONEI_CELL
RelaxSolverVolumeSteadyTest<tags::Cell, float> sse_relax_solver_vs_test_float_cell("VS CELL float");
RelaxSolverVolumeSteadyTest<tags::Cell, double> sse_relax_solver_vs_test_double_cell("VS CELL double");
#endif
