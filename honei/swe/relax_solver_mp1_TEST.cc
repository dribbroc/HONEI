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
#include <honei/swe/relax_solver.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>
#include <honei/swe/scenario.hh>
#include <honei/swe/scenario_manager.hh>
#include <sys/time.h>
using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;

template <typename Tag_, typename DT1_, typename DT2_>
class RelaxSolverMIXEDPRECTest :
    public BaseTest
{
    public:
        RelaxSolverMIXEDPRECTest(const std::string & type) :
            BaseTest("relax_solver_mixedprec_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long dwidth(40);
            unsigned long dheight(40);
            unsigned long timesteps(50);
            unsigned long k(15);

            DenseMatrix<DT1_> height(dheight, dwidth, DT1_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height.rows(); ++i)
            {
                for(ulint j=height.columns()-10; j<height.columns(); ++j)
                {
                    height[i][j] = DT1_(10);
                }
            }

            DenseMatrix<DT1_> bottom(dheight, dwidth, DT1_(1));

            DenseMatrix<DT1_> u1(dheight, dwidth, DT1_(0));
            DenseMatrix<DT1_> u2(dheight, dwidth, DT1_(0));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<DT1_> u(entries, DT1_(1));
            DenseVector<DT1_> v(entries, DT1_(1));
            DenseVector<DT1_> w(entries, DT1_(1)); 
            DenseVector<DT1_> bx (entries/3, DT1_(0));
            DenseVector<DT1_> by (entries/3, DT1_(0));
            DenseVector<DT1_> c (3,DT1_(5));
            DenseVector<DT1_> d (3,DT1_(5));
            c[0] = 10;
            c[1] = 6;
            c[2] = 11;
            d[0] = 10;
            d[1] = 5;
            d[2] = 11;

            DT1_ deltax = 5;
            DT1_ deltay = 5;
            DT1_ deltat = 5./24.;

            double eps = 10e-6;
            DT1_ manning = DT1_(0);
            ScenarioManager<DT1_, swe_solvers::RELAX, boundaries::REFLECT> scen_man;
            Scenario<DT1_, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);

            scen_man.allocate_scenario(&scenario);
            scen_man.allocate_scalarfields(&height, &bottom, &u1, &u2);
            scen_man.allocate_relax_vectors(&u, &v, &w, &c, &d);
            scen_man.allocate_bottom_slopes(&bx, &by);
            scen_man.set_environmental_variables(deltax, deltay, deltat, manning, eps);
//-------------------------------------------------------------------------------------------------------------------
            DenseMatrix<DT2_> height_2(dheight, dwidth, DT2_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height_2.rows(); ++i)
            {
                for(ulint j=height_2.columns()-10; j<height_2.columns(); ++j)
                {
                    height_2[i][j] = DT2_(10);
                }
            }

            DenseMatrix<DT2_> bottom_2(dheight, dwidth, DT2_(1));

            DenseMatrix<DT2_> u1_2(dheight, dwidth, DT2_(0));
            DenseMatrix<DT2_> u2_2(dheight, dwidth, DT2_(0));
            DenseVector<DT2_> u_2(entries, DT2_(1));
            DenseVector<DT2_> v_2(entries, DT2_(1));
            DenseVector<DT2_> w_2(entries, DT2_(1)); 
            DenseVector<DT2_> bx_2 (entries/3, DT2_(0));
            DenseVector<DT2_> by_2 (entries/3, DT2_(0));
            DenseVector<DT2_> c_2 (3,DT2_(5));
            DenseVector<DT2_> d_2 (3,DT2_(5));
            c_2[0] = 10;
            c_2[1] = 6;
            c_2[2] = 11;
            d_2[0] = 10;
            d_2[1] = 5;
            d_2[2] = 11;

            DT2_ deltax_2 = 5;
            DT2_ deltay_2 = 5;
            DT2_ deltat_2 = 5./24.;

            DT2_ manning_2 = DT2_(0);
            ScenarioManager<DT2_, swe_solvers::RELAX, boundaries::REFLECT> scen_man_2;
            Scenario<DT2_, swe_solvers::RELAX, boundaries::REFLECT> scenario_2(dwidth, dheight);

            scen_man_2.allocate_scenario(&scenario_2);
            scen_man_2.allocate_scalarfields(&height_2, &bottom_2, &u1_2, &u2_2);
            scen_man_2.allocate_relax_vectors(&u_2, &v_2, &w_2, &c_2, &d_2);
            scen_man_2.allocate_bottom_slopes(&bx_2, &by_2);
            scen_man_2.set_environmental_variables(deltax_2, deltay_2, deltat_2, manning_2, eps);

            ScenarioManager<DT2_, swe_solvers::RELAX, boundaries::REFLECT>::convert_scenario(scenario_2, scenario);

            if(scen_man_2.validate() && scen_man.validate())
            {
                RelaxSolver<Tag_, DT2_, DT2_, DT2_, DT2_, DT2_, source_types::SIMPLE, boundaries::REFLECT, precision_modes::FIXED> relax_solver (scenario_2);

                RelaxSolver<Tag_, DT1_, DT1_, DT1_, DT1_, DT1_, source_types::SIMPLE, boundaries::REFLECT, precision_modes::FIXED> relax_solver_double (scenario);
                relax_solver.do_preprocessing();
                string outHeight = stringify(height_2);

                for (ulint i = 1; i <= timesteps; ++i)
                {
                    if(i % k != 0)
                    {
                        relax_solver.solve();
                    }
                    else
                    {

                        ScenarioManager<DT1_, swe_solvers::RELAX, boundaries::REFLECT>::convert_scenario(scenario, scenario_2);
                        relax_solver_double.solve();
                        ScenarioManager<DT2_, swe_solvers::RELAX, boundaries::REFLECT>::convert_scenario(scenario_2, scenario);
                        //~relax_solver_double;
                    }
                    cout << "Timestep "<< i <<" / " << timesteps << " finished." <<endl;
                }

                bool pass = true;
                for(unsigned long i(0); i < height_2.rows(); ++i)
                {
                    for(unsigned long j(0); j < height_2.columns(); ++j)
                    {
                        if(height_2[i][j] < DT2_(4.8) || height_2[i][j] > DT2_(10.))
                            pass = false;
                    }
                }
                TEST_CHECK(pass);
            }
        }
};
RelaxSolverMIXEDPRECTest<tags::CPU, double, float> relax_solver_mp1_test("double to float");
RelaxSolverMIXEDPRECTest<tags::CPU::MultiCore, double, float> relax_solver_mp1_test_mc("MC double to float");
#ifdef HONEI_SSE
RelaxSolverMIXEDPRECTest<tags::CPU::SSE, double, float> relax_solver_mp1_test_sse("SSE double to float");
RelaxSolverMIXEDPRECTest<tags::CPU::MultiCore::SSE, double, float> relax_solver_mp1_test_mcsse("MC SSE double to float");
#endif
#ifdef HONEI_CELL
RelaxSolverMIXEDPRECTest<tags::Cell, double, float> relax_solver_mp1_test_cell("CELL double to float");
#endif
