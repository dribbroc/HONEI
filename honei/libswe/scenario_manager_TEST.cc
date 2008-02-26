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
#define SOLVER_BENCHMARK 1
#include <honei/libswe/relax_solver.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <string>
#include <honei/libswe/scenario.hh>
#include <honei/libswe/scenario_manager.hh>
#include <sys/time.h>

using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;

template <typename DataType_>
class ScenarioManagerTest :
    public BaseTest
{
    public:
        ScenarioManagerTest(const std::string & type) :
            BaseTest("scenario_manager_test<" + type + ">")
        {
            register_tag("CPU");
        }

        virtual void run() const
        {
            unsigned long dwidth = 40;
            unsigned long dheight = 40;
            unsigned long timesteps = 1;

            DenseMatrix<DataType_> height(dheight, dwidth, DataType_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height.rows(); ++i)
            {
                for(ulint j=height.columns()-10; j<height.columns(); ++j)
                {
                    height[i][j] = DataType_(10);
                }
            }

            DenseMatrix<DataType_> bottom(dheight, dwidth, DataType_(1));
            std::cout<<bottom<<std::endl;

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
            ScenarioManager<DataType_, swe_solvers::RELAX, boundaries::REFLECT> scen_man;
            Scenario<DataType_, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);

            scen_man.allocate_scenario(&scenario);
            scen_man.allocate_scalarfields(&height, &bottom, &u1, &u2);
            scen_man.allocate_relax_vectors(&u, &v, &w, &c, &d);
            scen_man.allocate_bottom_slopes(&bx, &by);
            scen_man.set_environmental_variables(deltax, deltay, deltat, manning, eps);

            if(scen_man.validate())
            {
                RelaxSolver<tags::CPU, DataType_, DataType_, DataType_, DataType_, DataType_, source_types::SIMPLE, boundaries::REFLECT, precision_modes::FIXED> relax_solver
                    (scenario);
                relax_solver.do_preprocessing();
                cout << "Height -field after preprocessing:\n";
                string outHeight = stringify(height);
                cout <<  outHeight;

                timeval start, end;
                for (ulint i = 1; i <= timesteps; ++i)
                {
                    gettimeofday(&start, 0);
                    relax_solver.solve();
                    gettimeofday(&end, 0);
                    cout << "Timestep "<< i <<" / " << timesteps << " finished." <<endl;
                    cout << "Solvetime: "<< end.tv_sec - start.tv_sec << " " << end.tv_usec - start.tv_usec<<endl;
                }
                cout << "Height -field after solve():\n";
                cout << stringify(height);

                bool pass = true;
                for(unsigned long i(0); i < height.rows(); ++i)
                {
                    for(unsigned long j(0); j < height.columns(); ++j)
                    {
                        if(height[i][j] < DataType_(4.8) || height[i][j] > DataType_(10.))
                            pass = false;
                    }
                }
                TEST_CHECK(pass);
            }
        }
};
ScenarioManagerTest<float> scenario_manager_test_float("float");

