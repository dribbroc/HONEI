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
            DenseMatrix<DataType_>* height = new DenseMatrix<DataType_> (3, 3, DataType_(1));
            DenseMatrix<DataType_>* bottom = new DenseMatrix<DataType_> (3, 3, DataType_(0));
            DenseMatrix<DataType_>* u1 = new DenseMatrix<DataType_> (3, 3, DataType_(1));
            DenseMatrix<DataType_>* u2 = new DenseMatrix<DataType_> (3, 3, DataType_(1));                                    
            DenseVector<DataType_>* u = new DenseVector<DataType_>(3*49, DataType_(0));
            DenseVector<DataType_>* v = new DenseVector<DataType_>(3*49, DataType_(0));
            DenseVector<DataType_>* w = new DenseVector<DataType_> (3*49, DataType_(0)); 
            DenseVector<DataType_> bx (49, DataType_(0));
            DenseVector<DataType_> by (49, DataType_(0));
            DenseVector<DataType_> c (3,DataType_(1));
            DenseVector<DataType_> d (3,DataType_(1));
            ulint dwith = 3;
            ulint dheight = 3;
            DataType_ deltax = 1;                       
            DataType_ deltay = 1;
            DataType_ deltat = 1;      
            double eps = 0.1;                  
            RelaxSolver<DataType_, DataType_, DataType_, DataType_, DataType_> relax_solver
                (height, bottom, u1, u2, u, v, w, 
                dwith, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d);
            relax_solver.do_preprocessing();
            cout << "Height -field after preprocessing:\n";
            string outHeight = stringify(*height);
            cout <<  outHeight;
            /*cout << "Relax - vectors after preprocessing:\n";
            cout << "u^T:\n";
            cout << stringify(*u) << endl;
            cout << "v^T:\n";
            cout << stringify(*v) << endl;
            cout << "w^T:\n";
            cout << stringify(*w) << endl;
            
            cout << "Bottom - slopes after preprocessing:\n";
            cout << "b_x^T:\n";
            cout << stringify(bx)<< endl;
            cout << "b_y^T:\n";
            cout << stringify(by)<< endl;
            */
            relax_solver.solve();
            cout << "Height -field after solve():\n";
            cout << stringify(*height);
            cout << "Relax - vectors after solve():\n";
            cout << "u^T:\n";
            cout << stringify(*u) << endl;
            cout << "v^T:\n";
            cout << stringify(*v) << endl;
            cout << "w^T:\n";
            cout << stringify(*w) << endl;
            delete height;
            delete bottom;
            delete u1;
            delete u2;
            delete u;
            delete v;
            delete w;
            TEST_CHECK(true);
        }
};
//RelaxSolverQuickTest<float> relax_solver_quick_test_float("float");
RelaxSolverQuickTest<double> relax_solver_quick_test_double("double");


