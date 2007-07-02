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
#include <libla/vector_scaled_sum.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_scaled_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(2)));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(3)));
                    
                DataType_ left(static_cast<DataType_>(2));
                DataType_ right(static_cast<DataType_>(3));
                                
                DenseVector<DataType_> sum1(VectorScaledSum<>::value(*dv1, *dv2, left, right));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 13 * size);
            }

            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(2,
                    static_cast<DataType_>(1)));
            DataType_ left(static_cast<DataType_>(2));
            DataType_ right(static_cast<DataType_>(3));

            TEST_CHECK_THROWS(VectorScaledSum<>::value(*dv00, *dv01, left, right), VectorSizeDoesNotMatch);
        }
};

DenseVectorScaledSumTest<float> dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<double> dense_vector_scaled__sum_test_double("double");

template <typename DataType_>
class DenseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_scaled_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                static_cast<DataType_>(2)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size,
                static_cast<DataType_>(3)));
                
            DataType_ left(static_cast<DataType_>(2));
            DataType_ right(static_cast<DataType_>(3));
                            
            DenseVector<DataType_> sum1(VectorScaledSum<>::value(*dv1, *dv2, left, right));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(v1, 13 * size);

            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(2,
                    static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorScaledSum<>::value(*dv00, *dv01, left, right), VectorSizeDoesNotMatch);
        }
};

DenseVectorScaledSumQuickTest<float> dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<double> dense_vector_scaled__sum_quick_test_double("double");

template <typename DataType_>
class SparseVectorScaledSumTest :
    public BaseTest
{
    public:
        SparseVectorScaledSumTest(const std::string & type) :
            BaseTest("sparse_vector_scaled_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > sv1(new SparseVector<DataType_>(size, size / 7 + 1));
                for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
                }            
                std::tr1::shared_ptr<SparseVector<DataType_> > sv2(new SparseVector<DataType_>(size, size / 8 + 1));
                for (typename Vector<DataType_>::ElementIterator i(sv2->begin_elements()), i_end(sv2->end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
                } 
                                     
                DataType_ left(static_cast<DataType_>(2));
                DataType_ right(static_cast<DataType_>(3));
                                
                SparseVector<DataType_> sum1(VectorScaledSum<>::value(*sv1, *sv2, left, right));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, size / 10 * 4 + size / 7 * 9);
            }

            std::tr1::shared_ptr<SparseVector<DataType_> > sv00(new SparseVector<DataType_>(1, 1));
            std::tr1::shared_ptr<SparseVector<DataType_> > sv01(new SparseVector<DataType_>(2, 1));
            DataType_ left00(static_cast<DataType_>(2));
            DataType_ right00(static_cast<DataType_>(3));

            TEST_CHECK_THROWS(VectorScaledSum<>::value(*sv00, *sv01, left00, right00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumTest<float> sparse_vector_scaled_sum_test_float("float");
SparseVectorScaledSumTest<double> sparse_vector_scaled__sum_test_double("double");

template <typename DataType_>
class SparseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_scaled_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(200);
            std::tr1::shared_ptr<SparseVector<DataType_> > sv1(new SparseVector<DataType_>(size, size / 6 + 1));
            for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
            }            
            std::tr1::shared_ptr<SparseVector<DataType_> > sv2(new SparseVector<DataType_>(size, size / 5 + 1));
            for (typename Vector<DataType_>::ElementIterator i(sv2->begin_elements()), i_end(sv2->end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
            } 
                                    
            DataType_ left(static_cast<DataType_>(2));
            DataType_ right(static_cast<DataType_>(3));
        
            std::tr1::shared_ptr<SparseVector<DataType_> > sum2(new SparseVector<DataType_>(size, size / 5 + 1));
            for (typename Vector<DataType_>::ElementIterator i(sum2->begin_elements()), i_end(sum2->end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(9);
                if (i.index() % 10 == 0) *i = static_cast<DataType_>(4);
                if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = static_cast<DataType_>(13);                
            }                          

            SparseVector<DataType_> sum1(VectorScaledSum<>::value(*sv1, *sv2, left, right));
//            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(sum1, *sum2);
//            TEST_CHECK_EQUAL(v1, (size / 10) * 4 + (size / 7) * 9);

            std::tr1::shared_ptr<SparseVector<DataType_> > sv00(new SparseVector<DataType_>(1, 1));
            std::tr1::shared_ptr<SparseVector<DataType_> > sv01(new SparseVector<DataType_>(2, 1));
            DataType_ left00(static_cast<DataType_>(2));
            DataType_ right00(static_cast<DataType_>(3));

            TEST_CHECK_THROWS(VectorScaledSum<>::value(*sv00, *sv01, left00, right00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumQuickTest<float> sparse_vector_scaled_sum_quick_test_float("float");
SparseVectorScaledSumQuickTest<double> sparse_vector_scaled__sum_quick_test_double("double");
