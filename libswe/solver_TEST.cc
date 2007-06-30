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

#include <string>

using namespace pg512;
using namespace tests;

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
            DenseMatrix<DataType_> height (3, 3, DataType_(1));
            DenseMatrix<DataType_> bottom (3, 3, DataType_(1));
            DenseMatrix<DataType_> u1 (3, 3, DataType_(1));
            DenseMatrix<DataType_> u2 (3, 3, DataType_(1));                                    
            DenseVector<DataType_> u (3*49, DataType_(1));
            DenseVector<DataType_> v(3*49, DataType_(1));
            DenseVector<DataType_> w (3*49, DataType_(1)); 
            DenseVector<DataType_> bx (49, DataType_(1));
            DenseVector<DataType_> by (49, DataType_(1));
            ulint dwith = 3;
            ulint dheight = 3;
            DataType_ deltax = 1;                       
            DataType_ deltay = 1;
            DataType_ deltat = 1;      
            double eps = 0.1;                  
            RelaxSolver<DataType_, DataType_, DataType_, DataType_, DataType_> relax_solver
                (&height, &bottom, &u1, &u2, &u, &v, &w, 
                dwith, dheight, deltax, deltay, deltat, eps, &bx, &by);
            relax_solver.do_preprocessing();
            TEST_CHECK(true);
        }
};
RelaxSolverQuickTest<float> relax_solver_quick_test_float("float");
RelaxSolverQuickTest<double> relax_solver_quick_test_double("double");


