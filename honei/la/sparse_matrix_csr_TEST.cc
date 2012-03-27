/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/la/sparse_matrix_csr.hh>
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
class SparseMatrixCSRQuickTest :
    public QuickTest
{
    public:
        SparseMatrixCSRQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_csr_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long blocks(1) ; blocks <= 4 ; blocks*=2)
            {
                Configuration::instance()->set_value("csr::blocksize", blocks);
                unsigned long size (111);
                SparseMatrix<DataType_> sms(size, size + 3);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sms.begin_elements()) ; i < sms.end_elements() ; ++i)
                {
                    if (i.index() % 5 == 0)
                        *i = (DataType_(i.index()) / 1.234 + 1);
                }
                sms[0][1] = 4711;
                sms[0][2] = 4711;
                sms[0][3] = 4711;

                SparseMatrixCSR<DataType_> sm0(sms);

                TEST_CHECK_EQUAL(sm0, sm0);
                TEST_CHECK_EQUAL(sm0, sm0.copy());
                TEST_CHECK_EQUAL(sm0.used_elements(), sms.used_elements());
                SparseMatrix<DataType_> sms2(sm0);
                TEST_CHECK_EQUAL(sms, sms2);
                SparseMatrixELL<DataType_> smell(sms);
                SparseMatrixELL<DataType_> smell2(sm0);
                TEST_CHECK_EQUAL(smell, smell2);
                SparseMatrixCSR<DataType_> sm1(smell);
                SparseMatrix<DataType_> sms3(sm1);
                TEST_CHECK_EQUAL(sms, sms3);
            }
        }
};
SparseMatrixCSRQuickTest<float> sparse_matrix_csr_quick_test_float("float");
SparseMatrixCSRQuickTest<double> sparse_matrix_csr_quick_test_double("double");

