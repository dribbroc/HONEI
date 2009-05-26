/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/vector_error.hh>
#include <unittest/unittest.hh>

#include <string>
#include <iostream>

using namespace honei;
using  namespace tests;

template <typename DataType_>
class SparseMatrixELLQuickTest :
    public QuickTest
{
    public:
        SparseMatrixELLQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_ell_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (5);
            DenseMatrix<DataType_> dm0(size, size + 3, DataType_(0));
            dm0(0, 0) = 1;
            dm0(0, 3) = 1;
            dm0(0, 4) = 1;
            dm0(2, 1) = 1;
            dm0(2, 3) = 1;
            dm0(3, 3) = 1;

            SparseMatrixELL<DataType_> sm0(dm0);

            std::cout<<dm0;
            std::cout<<sm0;

            TEST_CHECK_EQUAL(sm0, sm0);
            TEST_CHECK_EQUAL(sm0, sm0.copy());
        }
};
SparseMatrixELLQuickTest<float> sparse_matrix_ell_quick_test_float("float");
SparseMatrixELLQuickTest<double> sparse_matrix_ell_quick_test_double("double");

