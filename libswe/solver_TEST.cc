/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#include <libswe/solver.hh>
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <string>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class RelaxSolverQuickTest :
    public QuickTest
{
    public:
        RelaxSolverQuickTest(const std::string & type) :
            QuickTest("relax_solver_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            ulint dwidth =11;
            ulint dheight =11;
            ulint timesteps =1;

            DenseMatrix<DataType_> height(dheight, dwidth, DataType_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height.rows(); ++i)
            {
                for(ulint j=0; j<height.columns()-10; ++j)
                {
                    height[i][j] = DataType_(10);
                }
                 //(height)[0][i] = DataType_(10);
            }
            //END SCENARIO setup
            DenseMatrix<DataType_> bottom(dheight, dwidth, DataType_(1));
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
            DataType_ deltat = 5./22.;

            double eps = 10e-6;
            DataType_ manning = 0;

            RelaxSolver<Tag_, DataType_, DataType_, DataType_, DataType_, DataType_> relax_solver
                (&height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);
            relax_solver.do_preprocessing();
            cout << "Height -field after preprocessing:\n";
            string outHeight = stringify(height);
            cout <<  outHeight;
            for (ulint i = 1; i <= timesteps; ++i)
            {
                relax_solver.solve();
            }
            cout << "Height -field after solve():\n";
            cout << stringify(height);
            /*cout << "Relax - vectors after solve():\n";
            cout << "u^T:\n";
            cout << stringify(u) << endl;
            cout << "v^T:\n";
            cout << stringify(v) << endl;
            cout << "w^T:\n";
            cout << stringify(w) << endl;*/
            /*delete height;
            delete bottom;
            delete u1;
            delete u2;
            delete u;
            delete v;
            delete w;*/
            TEST_CHECK(true);
        }
};
//RelaxSolverQuickTest<float> relax_solver_quick_test_float("float");
//#ifndef HONEI_SSE
RelaxSolverQuickTest<tags::CPU, double> relax_solver_quick_test_double("double");
/*#endif
#ifdef HONEI_SSE
RelaxSolverQuickTest<tags::CPU::SSE, double> relax_solver_quick_test_double("sse double");
#endif*/
