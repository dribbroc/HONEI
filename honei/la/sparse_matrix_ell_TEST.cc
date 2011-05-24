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
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/unittest.hh>
#include <honei/util/configuration.hh>

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
            for (unsigned long threads(1) ; threads <= 16 ; threads*=2)
            {
                Configuration::instance()->set_value("ell::threads", threads);
                unsigned long size (111);
                SparseMatrix<DataType_> sms(size, size + 3);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sms.begin_elements()) ; i < sms.end_elements() ; ++i)
                {
                    if (i.index() % 5 == 0)
                        *i = (DataType_(i.index()) / 1.234 + 1);
                }
                sms[0][0] = 4711;

                SparseMatrixELL<DataType_> sm0(sms);

                TEST_CHECK_EQUAL(sm0, sm0);
                TEST_CHECK_EQUAL(sm0, sm0.copy());
                TEST_CHECK_EQUAL(sms.used_elements(), sm0.used_elements());
                for (unsigned long row(0) ; row < sm0.rows() ; ++row)
                {
                    for (unsigned long col(0) ; col < sm0.columns() ; ++col)
                    {
                        TEST_CHECK_EQUAL(sm0(row, col), sms(row, col));
                    }
                }
                SparseMatrix<DataType_> sms2(sm0);
                TEST_CHECK_EQUAL(sms, sms2);
            }
        }
};
SparseMatrixELLQuickTest<float> sparse_matrix_ell_quick_test_float("float");
SparseMatrixELLQuickTest<double> sparse_matrix_ell_quick_test_double("double");

