/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <iostream>

using namespace pg512;
using namespace tests;
using namespace std;

template <typename DataType_>
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
            ulint dwidth =2;
            ulint dheight =2;
            ulint timesteps = 1;
 
            DenseMatrix<DataType_>* height = new DenseMatrix<DataType_> (dheight, dwidth, DataType_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height->rows(); ++i)
            {
                (*height)[0][i] = 5;
            }
            //END SCENARIO setup
            DenseMatrix<DataType_>* bottom = new DenseMatrix<DataType_> (dheight, dwidth, DataType_(1));
            DenseMatrix<DataType_>* u1 = new DenseMatrix<DataType_> (dheight, dwidth, DataType_(1));
            DenseMatrix<DataType_>* u2 = new DenseMatrix<DataType_> (dheight, dwidth, DataType_(1));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<DataType_>* u = new DenseVector<DataType_>(entries, DataType_(1));
            DenseVector<DataType_>* v = new DenseVector<DataType_>(entries, DataType_(1));
            DenseVector<DataType_>* w = new DenseVector<DataType_>(entries, DataType_(1)); 
            DenseVector<DataType_> bx (entries/3, DataType_(0));
            DenseVector<DataType_> by (entries/3, DataType_(0));
            DenseVector<DataType_> c (3,DataType_(1));
            DenseVector<DataType_> d (3,DataType_(1));
            
            DataType_ deltax = 1;                       
            DataType_ deltay = 1;
            DataType_ deltat = 0.5;

            double eps = 0.1;                  
            RelaxSolver<DataType_, DataType_, DataType_, DataType_, DataType_> relax_solver
                (height, bottom, u1, u2, u, v, w, 
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d);
            relax_solver.do_preprocessing();
            cout << "Height -field after preprocessing:\n";
            string outHeight = stringify(*height);
            cout <<  outHeight;
            for (ulint i =1; i <= timesteps; ++i) 
            {
                relax_solver.solve();
            }
            cout << "Height -field after solve():\n";
            cout << stringify(*height);
            /*cout << "Relax - vectors after solve():\n";
            cout << "u^T:\n";
            cout << stringify(*u) << endl;
            cout << "v^T:\n";
            cout << stringify(*v) << endl;
            cout << "w^T:\n";
            cout << stringify(*w) << endl;*/
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
RelaxSolverQuickTest<double> relax_solver_quick_test_double("double");


