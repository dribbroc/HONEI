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
#include <libla/vector_sum.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));

                DenseVector<DataType_> sum1(VectorSum<>::value(*dv1, *dv2));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 2 * size);
            }

            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(2,
                    static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorSum<>::value(*dv00, *dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorSumTest<float> dense_vector_sum_test_float("float");
DenseVectorSumTest<double> dense_vector_sum_test_double("double");

template <typename DataType_>
class DenseVectorSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size,
                static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                static_cast<DataType_>(1)));

            DenseVector<DataType_> sum1(VectorSum<>::value(*dv1, *dv2));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(v1, 2 * size);

            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(2,
                    static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorSum<>::value(*dv00, *dv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorSumQuickTest<float>  dense_vector_sum_quick_test_float("float");
DenseVectorSumQuickTest<double> dense_vector_sum_quick_test_double("double");

template <typename DataType_>
class SparseVectorSumTest :
    public BaseTest
{
    public:
        SparseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > sv1(new SparseVector<DataType_>(size, size / 5 + 1)),
                    sv2(new SparseVector<DataType_>(size, size / 5 + 1));
                for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()),
                    j(sv2->begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = static_cast<DataType_>(1);
                        *j = static_cast<DataType_>(2);
                    }
                }    
                SparseVector<DataType_> sum1(VectorSum<>::value(*sv1, *sv2));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 3 * size / 10);
            }

            std::tr1::shared_ptr<SparseVector<DataType_> > sv00(new SparseVector<DataType_>(1, 1));
            std::tr1::shared_ptr<SparseVector<DataType_> > sv01(new SparseVector<DataType_>(2, 1));

            TEST_CHECK_THROWS(VectorSum<>::value(*sv00, *sv01), VectorSizeDoesNotMatch);
        }
};

SparseVectorSumTest<float> sparse_vector_sum_test_float("float");
SparseVectorSumTest<double> sparse_vector_sum_test_double("double");

template <typename DataType_>
class SparseVectorSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            std::tr1::shared_ptr<SparseVector<DataType_> > sv1(new SparseVector<DataType_>(size, size / 5 + 1)),
                sv2(new SparseVector<DataType_>(size, size / 5 + 1));
            for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()),
                j(sv2->begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = static_cast<DataType_>(1);
                    *j = static_cast<DataType_>(2);
                }
            }    
            SparseVector<DataType_> sum1(VectorSum<>::value(*sv1, *sv2));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(v1, 3 * size / 10);

            std::tr1::shared_ptr<SparseVector<DataType_> > sv00(new SparseVector<DataType_>(1, 1));
            std::tr1::shared_ptr<SparseVector<DataType_> > sv01(new SparseVector<DataType_>(2, 1));

            TEST_CHECK_THROWS(VectorSum<>::value(*sv00, *sv01), VectorSizeDoesNotMatch);
        }
};

SparseVectorSumQuickTest<float> sparse_vector_sum_quick_test_float("float");
SparseVectorSumQuickTest<double> sparse_vector_sum_quick_test_double("double");
