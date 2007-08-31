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

#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_masked_min_index.hh>
#include <libla/vector_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using  namespace tests;

template <typename DataType_>
class DenseVectorMaskedMinIndexTest :
    public BaseTest
{
    public:
        DenseVectorMaskedMinIndexTest(const std::string & type) :
            BaseTest("dense_vector_masked_min_index_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                DenseVector<bool> mask(size, false);
                DenseVector<DataType_> dv(size);

                for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()),
                    i_end(dv.end_elements()) ; i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                (dv)[2] = DataType_(0.009);

                for (typename Vector<bool>::ElementIterator i(mask.begin_elements()),
                    i_end(mask.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 2 == 0)
                    {
                        *i=true;
                    }
                }
                unsigned long result(VectorMaskedMinIndex<DataType_>::value(dv, mask));
                TEST_CHECK_EQUAL(result, 2);
                (dv)[2] = DataType_(0);
                unsigned long result2(VectorMaskedMinIndex<DataType_>::value(dv, mask));
                TEST_CHECK_EQUAL(result2, 2);
                (dv)[2] = DataType_(-5.2345);
                unsigned long result3(VectorMaskedMinIndex<DataType_>::value(dv, mask));
                TEST_CHECK_EQUAL(result3, 2);
            }

            DenseVector<bool> mask1(2, false);
            DenseVector<DataType_> dv1(3, DataType_(1));

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(dv1, mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMinIndexTest<float> dense_vector_masked_min_index_test_float("float");
DenseVectorMaskedMinIndexTest<double> dense_vector_masked_min_index_test_double("double");

template <typename DataType_>
class DenseVectorMaskedMinIndexQuickTest :
    public QuickTest
{
    public:
        DenseVectorMaskedMinIndexQuickTest(const std::string & type) :
            QuickTest("dense_vector_masked_min_index_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DataType_ count = 0;
            DenseVector<bool> mask(size, false);
            DenseVector<DataType_> dv(size);

            for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()),
                i_end(dv.end_elements()) ; i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }

            (dv)[2] = DataType_(0.009);

            for (typename Vector<bool>::ElementIterator i(mask.begin_elements()),
                i_end(mask.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 2 == 0)
                {
                    *i=true;
                }
            }
            unsigned long result(VectorMaskedMinIndex<DataType_>::value(dv, mask));
            TEST_CHECK_EQUAL(result, 2);
            (dv)[2] = DataType_(0);
            unsigned long result2(VectorMaskedMinIndex<DataType_>::value(dv, mask));
            TEST_CHECK_EQUAL(result2, 2);
            (dv)[2] = DataType_(-5.2345);
            unsigned long result3(VectorMaskedMinIndex<DataType_>::value(dv, mask));
            TEST_CHECK_EQUAL(result3, 2);


            DenseVector<bool> mask1(2, false);
            DenseVector<DataType_> dv1(3, DataType_(1));

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(dv1, mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMinIndexQuickTest<float>  dense_vector_masked_min_index_quick_test_float("float");
DenseVectorMaskedMinIndexQuickTest<double> dense_vector_masked_min_index_quick_test_double("double");

template <typename DataType_>
class SparseVectorMaskedMinIndexTest :
    public BaseTest
{
    public:
        SparseVectorMaskedMinIndexTest(const std::string & type) :
            BaseTest("sparse_vector_masked_min_index_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                DenseVector<bool> mask(size, false);

                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }                 

                for (typename Vector<bool>::ElementIterator i(mask.begin_elements()),
                    i_end(mask.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 2 == 0)
                    {
                        *i=true;
                    }
                }
                (sv1)[2] = DataType_(0);
                unsigned long result2(VectorMaskedMinIndex<DataType_>::value(sv1, mask));
                TEST_CHECK_EQUAL(result2, 2);
                (sv1)[2] = DataType_(-5.2345);
                unsigned long result3(VectorMaskedMinIndex<DataType_>::value(sv1, mask));
                TEST_CHECK_EQUAL(result3, 2);
            }

            DenseVector<bool> mask01(2, false);
            SparseVector<DataType_> sv01(3, 1);

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(sv01, mask01), VectorSizeDoesNotMatch);
        }
};
SparseVectorMaskedMinIndexTest<float> sparse_vector_masked_min_index_test_float("float");
SparseVectorMaskedMinIndexTest<double> sparse_vector_masked_min_index_test_double("double");

template <typename DataType_>
class SparseVectorMaskedMinIndexQuickTest :
    public QuickTest
{
    public:
        SparseVectorMaskedMinIndexQuickTest(const std::string & type) :
            QuickTest("sparse_vector_masked_min_index_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DataType_ count = 0;
            DenseVector<bool> mask(size, false);

            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }                 

            for (typename Vector<bool>::ElementIterator i(mask.begin_elements()),
                i_end(mask.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 2 == 0)
                {
                    *i=true;
                }
            }
            (sv1)[2] = DataType_(0);
            unsigned long result2(VectorMaskedMinIndex<DataType_>::value(sv1, mask));
            TEST_CHECK_EQUAL(result2, 2);
            (sv1)[2] = DataType_(-5.2345);
            unsigned long result3(VectorMaskedMinIndex<DataType_>::value(sv1, mask));
            TEST_CHECK_EQUAL(result3, 2);


            DenseVector<bool> mask01(2, false);
            SparseVector<DataType_> sv01(3, 1);

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(sv01, mask01), VectorSizeDoesNotMatch);
        }
};
SparseVectorMaskedMinIndexQuickTest<float> sparse_vector_masked_min_index_quick_test_float("float");
SparseVectorMaskedMinIndexQuickTest<double> sparse_vector_masked_min_index_quick_test_double("double");
