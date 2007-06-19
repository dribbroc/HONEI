/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
 
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_row_sum_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseMatrixRowSumVectorTest :
    public BaseTest
{
    public:
        DenseMatrixRowSumVectorTest(const std::string & type) :
            BaseTest("dense_matrix_row_sum_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index() + 1);
                }

                std::tr1::shared_ptr<DenseVector<DataType_> > dv(MatrixRowSumVector<DataType_>::value(*dm));
                DataType_ s(size);
                for (typename Vector<DataType_>::ElementIterator v(dv->begin_elements()), v_end(dv->end_elements()) ;
                        v != v_end ; ++v)
                {
                    DataType_ last((v.index() + 1) * s), first(v.index() * s);
                    DataType_ ssum((last * (last + 1) / 2) - (first * (first + 1) / 2));

                    TEST_CHECK_EQUAL_WITHIN_EPS(*v, ssum, ssum * 200 * std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseMatrixRowSumVectorTest<float> dense_matrix_row_sum_vector_test_float("float");
DenseMatrixRowSumVectorTest<double> dense_matrix_row_sum_vector_test_double("double");



template <typename DataType_>
class DenseMatrixRowSumVectorQuickTest :
    public QuickTest
{
    public:
        DenseMatrixRowSumVectorQuickTest(const std::string & type) :
            QuickTest("dense_matrix_row_sum_vector_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(3);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size+2));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>(i.index() + 1);
            }

            std::tr1::shared_ptr<DenseVector<DataType_> > dv(MatrixRowSumVector<DataType_>::value(*dm));
            DataType_ s(size);
            for (typename Vector<DataType_>::ElementIterator v(dv->begin_elements()), v_end(dv->end_elements()) ;
                    v != v_end ; ++v)
            {
                DataType_ last((v.index() + 1) * s), first(v.index() * s);
                DataType_ ssum((last * (last + 1) / 2) - (first * (first + 1) / 2));

                TEST_CHECK_EQUAL_WITHIN_EPS(*v, ssum, ssum * 100 * std::numeric_limits<DataType_>::epsilon());
            }
        }
};
DenseMatrixRowSumVectorQuickTest<float>  dense_matrix_row_sum_vector_quick_test_float("float");
DenseMatrixRowSumVectorQuickTest<double> dense_matrix_row_sum_vector_quick_test_double("double");
