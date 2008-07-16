/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Joachim Messer <joachim.messer@uni-dortmund.de>
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

#include <honei/la/dense_vector.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/element_product.hh>
#include <honei/la/vector_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class DenseVectorElementProductTest :
    public BaseTest
{
    public:
        DenseVectorElementProductTest(const std::string & type) :
            BaseTest("dense_vector_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0)),
                    dv3(size, DataType_(0));

                unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2
                DataType_ sign_1(1), sign_2(1);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    *i = sign_1 * (i.index() % num_limit);
                    *j = sign_2 * (num_limit - (i.index() % num_limit) - 1);
                    *k = sign_1 * sign_2 * ((i.index() % num_limit) * (num_limit - (i.index() % num_limit) - 1));
                    sign_1 *= -1;
                    sign_2 *= (-1) * sign_1;
                }
                ElementProduct<Tag_>::value(dv1, dv2);

                TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
            }

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dv02, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorElementProductTest<tags::CPU, float> dense_vector_elementwise_product_test_float("float");
DenseVectorElementProductTest<tags::CPU, double> dense_vector_elementwise_product_test_double("double");
DenseVectorElementProductTest<tags::CPU::MultiCore, float> mc_dense_vector_elementwise_product_test_float("MC float");
DenseVectorElementProductTest<tags::CPU::MultiCore, double> mc_dense_vector_elementwise_product_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorElementProductTest<tags::CPU::SSE, float> sse_dense_vector_elementwise_product_test_float("sse float");
DenseVectorElementProductTest<tags::CPU::SSE, double> sse_dense_vector_elementwise_product_test_double("sse double");
DenseVectorElementProductTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_elementwise_product_test_float("MC SSE float");
DenseVectorElementProductTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_elementwise_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorElementProductTest<tags::GPU::CUDA, float> cuda_dense_vector_elementwise_product_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorElementProductTest<tags::Cell, float> cell_dense_vector_element_product_test_float("Cell float");
DenseVectorElementProductTest<tags::Cell, double> cell_dense_vector_element_product_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class DenseVectorElementProductQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(4711);
            DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0)),
                dv3(size, DataType_(0));

            unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2
            DataType_ sign_1(1), sign_2(1);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                j(dv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                *i = sign_1 * (i.index() % num_limit);
                *j = sign_2 * (num_limit - (i.index() % num_limit) - 1);
                *k = sign_1 * sign_2 * ((i.index() % num_limit) * (num_limit - (i.index() % num_limit) - 1));
                sign_1 *= -1;
                sign_2 *= (-1) * sign_1;
            }

            ElementProduct<Tag_>::value(dv1, dv2);

            TEST(dv1.lock(lm_read_only),
                    for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        k(dv3.begin_elements()) ; i != i_end ; ++i, ++k)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *k, std::numeric_limits<DataType_>::epsilon());
                    },
                    dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dv02, dv01), VectorSizeDoesNotMatch);

        }
};
DenseVectorElementProductQuickTest<tags::CPU, float> dense_vector_elementwise_product_quick_test_float("float");
DenseVectorElementProductQuickTest<tags::CPU, double> dense_vector_elementwise_product_quick_test_double("double");
DenseVectorElementProductQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_elementwise_product_quick_test_float("MC float");
DenseVectorElementProductQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_elementwise_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorElementProductQuickTest<tags::CPU::SSE, float> sse_dense_vector_elementwise_product_quick_test_float("sse float");
DenseVectorElementProductQuickTest<tags::CPU::SSE, double> sse_dense_vector_elementwise_product_quick_test_double("sse double");
DenseVectorElementProductQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_elementwise_product_quick_test_float("MC SSE float");
DenseVectorElementProductQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_elementwise_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorElementProductQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_elementwise_product_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorElementProductQuickTest<tags::Cell, float> cell_dense_vector_element_product_quick_test_float("Cell float");
DenseVectorElementProductQuickTest<tags::Cell, double> cell_dense_vector_element_product_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeElementProductTest :
    public BaseTest
{
    public:
        DenseVectorRangeElementProductTest(const std::string & type) :
            BaseTest("dense_vector_range_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0)),
                    dv3(size, DataType_(0));

                unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2
                DataType_ sign_1(1), sign_2(1);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()); i != i_end ; ++i, ++j)
                {
                    *i = sign_1 * (i.index() % num_limit);
                    *j = sign_2 * (num_limit - (i.index() % num_limit) - 1);
                    sign_1 *= -1;
                    sign_2 *= (-1) * sign_1;
                }

                for (unsigned long range_size(1); range_size <= size-8; ++range_size)
                {
                    for (unsigned long offset1(4); offset1 < 8; ++offset1)
                    {
                        for (unsigned long offset2(4); offset2 < 8; ++offset2)
                        {
                            DenseVectorRange<DataType_> dvr1(dv1.copy(), range_size, offset1);
                            DenseVectorRange<DataType_> dvr2(dv2, range_size, offset2);
                            DenseVector<DataType_> dvr3(range_size, DataType_(0));
                            for (typename Vector<DataType_>::ElementIterator i(dvr3.begin_elements()), 
                                i_end(dvr3.end_elements()); i != i_end ; ++i)
                            {
                                *i = ((i.index() + offset1)%2? DataType_(-1) : DataType_(1)) * ((i.index() + offset2)%4 < 2? DataType_(1) : DataType_(-1)) * ((i.index()+offset1)%num_limit) * (num_limit - ((i.index() + offset2)%num_limit) -1);
                            }

                            dvr1.lock(lm_read_and_write);
                            dvr1.unlock(lm_read_and_write);
                            dvr2.lock(lm_read_and_write);
                            dvr2.unlock(lm_read_and_write);
                            ElementProduct<Tag_>::value(dvr1, dvr2);
                            TEST(dvr1.lock(lm_read_only), TEST_CHECK_EQUAL(dvr1, dvr3), dvr1.unlock(lm_read_only));
                        }
                    }
                }

            }

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dv02, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorRangeElementProductTest<tags::CPU, float> dense_vector_range_elementwise_product_test_float("float");
DenseVectorRangeElementProductTest<tags::CPU, double> dense_vector_range_elementwise_product_test_double("double");
DenseVectorRangeElementProductTest<tags::CPU::MultiCore, float> mc_dense_vector_range_elementwise_product_test_float("MC float");
DenseVectorRangeElementProductTest<tags::CPU::MultiCore, double> mc_dense_vector_range_elementwise_product_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeElementProductTest<tags::CPU::SSE, float> sse_dense_vector_range_elementwise_product_test_float("sse float");
DenseVectorRangeElementProductTest<tags::CPU::SSE, double> sse_dense_vector_range_elementwise_product_test_double("sse double");
DenseVectorRangeElementProductTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_elementwise_product_test_float("MC SSE float");
DenseVectorRangeElementProductTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_elementwise_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeElementProductTest<tags::GPU::CUDA, float> cuda_dense_vector_range_elementwise_product_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeElementProductTest<tags::Cell, float> cell_dense_vector_range_element_product_test_float("Cell float");
DenseVectorRangeElementProductTest<tags::Cell, double> cell_dense_vector_range_element_product_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeElementProductQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeElementProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
#if defined HONEI_CELL
            unsigned long size(10000);
#else
            unsigned long size(471);
#endif
            DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0)),
                dv3(size, DataType_(0));

            unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2
            DataType_ sign_1(1), sign_2(1);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                j(dv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                *i = sign_1 * (i.index() % num_limit);
                *j = sign_2 * (num_limit - (i.index() % num_limit) - 1);
                sign_1 *= -1;
                sign_2 *= (-1) * sign_1;
            }
            for (unsigned long range_size(size-16); range_size <= size-8; ++range_size)
            {
                for (unsigned long offset1(4); offset1 < 8; ++offset1)
                {
                    for (unsigned long offset2(4); offset2 < 8; ++offset2)
                    {
                        DenseVectorRange<DataType_> dvr1(dv1.copy(), range_size, offset1);
                        DenseVectorRange<DataType_> dvr2(dv2, range_size, offset2);
                        DenseVector<DataType_> dvr3(range_size, DataType_(0));
                        for (typename Vector<DataType_>::ElementIterator i(dvr3.begin_elements()),
                            i_end(dvr3.end_elements()); i != i_end ; ++i)
                        {
                            *i = ((i.index() + offset1)%2? DataType_(-1) : DataType_(1)) * ((i.index() + offset2)%4 < 2? DataType_(1) : DataType_(-1)) * ((i.index()+offset1)%num_limit) * (num_limit - ((i.index()+offset2)%num_limit) -1);
                        }

                        dvr1.lock(lm_read_and_write);
                        dvr1.unlock(lm_read_and_write);
                        dvr2.lock(lm_read_and_write);
                        dvr2.unlock(lm_read_and_write);
                        ElementProduct<Tag_>::value(dvr1, dvr2);
                        TEST(dvr1.lock(lm_read_only), TEST_CHECK_EQUAL(dvr1, dvr3), dvr1.unlock(lm_read_only));
                    }
                }
            }

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dv02, dv01), VectorSizeDoesNotMatch);

        }
};
DenseVectorRangeElementProductQuickTest<tags::CPU, float> dense_vector_range_elementwise_product_quick_test_float("float");
DenseVectorRangeElementProductQuickTest<tags::CPU, double> dense_vector_range_elementwise_product_quick_test_double("double");
DenseVectorRangeElementProductQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_range_elementwise_product_quick_test_float("MC float");
DenseVectorRangeElementProductQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_elementwise_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeElementProductQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_elementwise_product_quick_test_float("sse float");
DenseVectorRangeElementProductQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_elementwise_product_quick_test_double("sse double");
DenseVectorRangeElementProductQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_elementwise_product_quick_test_float("MC SSE float");
DenseVectorRangeElementProductQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_elementwise_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeElementProductQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_range_elementwise_product_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeElementProductQuickTest<tags::Cell, float> cell_dense_vector_range_element_product_quick_test_float("Cell float");
DenseVectorRangeElementProductQuickTest<tags::Cell, double> cell_dense_vector_range_element_product_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class SparseVectorDenseVectorElementProductTest :
    public BaseTest
{
    public:
        SparseVectorDenseVectorElementProductTest(const std::string & type) :
            BaseTest("sparse_vector_dense_vector_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = static_cast<DataType_>(2);
                }
                DenseVector<DataType_> dv2(size, DataType_(3));
                SparseVector<DataType_> sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0)
                        *i = static_cast<DataType_>(6);
                }

                ElementProduct<Tag_>::value(sv1, dv2);

                TEST_CHECK_EQUAL(sv1, sv3);
            }

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementProductTest<tags::CPU, float> sparse_vector_dense_vector_elementwise_product_test_float("float");
SparseVectorDenseVectorElementProductTest<tags::CPU, double> sparse_vector_dense_vector_elementwise_product_test_double("double");
SparseVectorDenseVectorElementProductTest<tags::CPU::MultiCore, float> mc_sparse_vector_dense_vector_elementwise_product_test_float("MC float");
SparseVectorDenseVectorElementProductTest<tags::CPU::MultiCore, double> mc_sparse_vector_dense_vector_elementwise_product_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseVectorDenseVectorElementProductQuickTest :
    public QuickTest
{
    public:
        SparseVectorDenseVectorElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_vector_dense_vector_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(2);
            }
            DenseVector<DataType_> dv2(size, DataType_(3));
            SparseVector<DataType_> sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0)
                    *i = static_cast<DataType_>(6);
            }
            ElementProduct<Tag_>::value(sv1, dv2);

            TEST_CHECK_EQUAL(sv1, sv3);

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementProductQuickTest<tags::CPU, float> sparse_vector_dense_vector_elementwise_product_quick_test_float("float");
SparseVectorDenseVectorElementProductQuickTest<tags::CPU, double> sparse_vector_dense_vector_elementwise_product_quick_test_double("double");
SparseVectorDenseVectorElementProductQuickTest<tags::CPU::MultiCore, float> mc_sparse_vector_dense_vector_elementwise_product_quick_test_float("MC float");
SparseVectorDenseVectorElementProductQuickTest<tags::CPU::MultiCore, double> mc_sparse_vector_dense_vector_elementwise_product_quick_test_double("MC double");


template <typename DataType_>
class SparseVectorElementProductTest :
    public BaseTest
{
    public:
        SparseVectorElementProductTest(const std::string & type) :
            BaseTest("sparse_vector_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
                }
                SparseVector<DataType_> sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0 && i.index() % 7 == 0)
                        *i = static_cast<DataType_>(6);
                }
                ElementProduct<>::value(sv1, sv2);

                TEST_CHECK_EQUAL(sv1, sv3);
            }

            SparseVector<DataType_> sv01(3, 1), sv02(4, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sv02, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorElementProductTest<float> sparse_vector_elementwise_product_test_float("float");
SparseVectorElementProductTest<double> sparse_vector_elementwise_product_test_double("double");

template <typename DataType_>
class SparseVectorElementProductQuickTest :
    public QuickTest
{
    public:
        SparseVectorElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_vector_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(100);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
            }
            SparseVector<DataType_> sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                    *i = static_cast<DataType_>(6);
            }
            ElementProduct<>::value(sv1, sv2);

            TEST_CHECK_EQUAL(sv1, sv3);

            SparseVector<DataType_> sv01(3, 1), sv02(4, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sv02, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorElementProductQuickTest<float> sparse_vector_elementwise_product_quick_test_float("float");
SparseVectorElementProductQuickTest<double> sparse_vector_elementwise_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class BandedMatrixDenseMatrixElementProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixElementProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size, DataType_(3));
                DenseVector<DataType_> dv2(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3);

                ElementProduct<Tag_>::value(bm2, dm1);

                TEST_CHECK_EQUAL(bm2, bm3);
            }

            DenseMatrix<DataType_> dm01(6, 5), dm03(5, 6);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, dm01), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixDenseMatrixElementProductTest<tags::CPU, float> banded_matrix_dense_matrix_elementwise_product_test_float("float");
BandedMatrixDenseMatrixElementProductTest<tags::CPU, double> banded_matrix_dense_matrix_elementwise_product_test_double("double");

template <typename Tag_, typename DataType_>
class BandedMatrixDenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size, DataType_(3));
            DenseVector<DataType_> dv2(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3);

            ElementProduct<Tag_>::value(bm2, dm1);

            TEST_CHECK_EQUAL(bm2, bm3);

            DenseMatrix<DataType_> dm01(6, 5), dm03(5, 6);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, dm01), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixDenseMatrixElementProductQuickTest<tags::CPU, float> banded_matrix_dense_matrix_elementwise_product_quick_test_float("float");
BandedMatrixDenseMatrixElementProductQuickTest<tags::CPU, double> banded_matrix_dense_matrix_elementwise_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class BandedMatrixElementProductTest :
    public BaseTest
{
    public:
        BandedMatrixElementProductTest(const std::string & type) :
            BaseTest("banded_matrix_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));
                DenseVector<DataType_> dv2(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);

                ElementProduct<Tag_>::value(bm1, bm2);

                TEST_CHECK_EQUAL(bm1, bm3);
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixElementProductTest<tags::CPU, float> banded_matrix_elementwise_product_test_float("float");
BandedMatrixElementProductTest<tags::CPU, double> banded_matrix_elementwise_product_test_double("double");
BandedMatrixElementProductTest<tags::CPU::MultiCore, float> mc_banded_matrix_elementwise_product_test_float("MC float");
BandedMatrixElementProductTest<tags::CPU::MultiCore, double> mc_banded_matrix_elementwise_product_test_double("MC double");

template <typename Tag_, typename DataType_>
class BandedMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> dv1(size, DataType_(3));
            DenseVector<DataType_> dv2(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);

            ElementProduct<Tag_>::value(bm1, bm2);

            TEST_CHECK_EQUAL(bm1, bm3);

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixElementProductQuickTest<tags::CPU, float> banded_matrix_elementwise_product_quick_test_float("float");
BandedMatrixElementProductQuickTest<tags::CPU, double> banded_matrix_elementwise_product_quick_test_double("double");
BandedMatrixElementProductQuickTest<tags::CPU::MultiCore, float> mc_banded_matrix_elementwise_product_quick_test_float("MC float");
BandedMatrixElementProductQuickTest<tags::CPU::MultiCore, double> mc_banded_matrix_elementwise_product_quick_test_double("MC double");

template <typename DataType_>
class BandedMatrixSparseMatrixElementProductTest :
    public BaseTest
{
    public:
        BandedMatrixSparseMatrixElementProductTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));
                BandedMatrix<DataType_> bm1(size, dv1),  bm3(size);
                SparseMatrix<DataType_> sm2(size, size, size / 7 + 1);
                typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements());
                typename BandedMatrix<DataType_>::ElementIterator k(bm3.begin_elements());
                for (typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements()),
                    j_end(sm2.end_elements()) ; j != j_end ; ++i, ++j, ++k)
                {
                    if (j.index() % 7 == 0)
                        *j = DataType_(2);

                    if (*i != DataType_(0) && i.index() % 7 == 0)
                        *k = DataType_(6);
                }

                ElementProduct<>::value(bm1, sm2);

                TEST_CHECK_EQUAL(bm1, bm3);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm02), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm03), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixSparseMatrixElementProductTest<float> banded_matrix_sparse_matrix_elementwise_product_test_float("float");
BandedMatrixSparseMatrixElementProductTest<double> banded_matrix_sparse_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(10);
            DenseVector<DataType_> dv1(size, DataType_(3));
            BandedMatrix<DataType_> bm1(size, dv1),  bm3(size);
            SparseMatrix<DataType_> sm2(size, size, size / 7 + 1);
            typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements());
            typename BandedMatrix<DataType_>::ElementIterator k(bm3.begin_elements());
            for (typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements()),
                j_end(sm2.end_elements()) ; j != j_end ; ++i, ++j, ++k)
            {
                if (j.index() % 7 == 0)
                    *j = DataType_(2);

                if (*i != DataType_(0) && i.index() % 7 == 0)
                    *k = DataType_(6);
            }

            ElementProduct<>::value(bm1, sm2);

            TEST_CHECK_EQUAL(bm1, bm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm02), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm03), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixSparseMatrixElementProductQuickTest<float> banded_matrix_sparse_matrix_elementwise_product_quick_test_float("float");
BandedMatrixSparseMatrixElementProductQuickTest<double> banded_matrix_sparse_matrix_elementwise_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class DenseMatrixElementProductTest :
    public BaseTest
{
    public:
        DenseMatrixElementProductTest(const std::string & type) :
            BaseTest("dense_matrix_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(0)), dm2(size+1, size, DataType_(0)),
                    dm3(size+1, size, DataType_(0));

                unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2
                DataType_ sign_1(1), sign_2(1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), 
                    i_end(dm1.end_elements()), j(dm2.begin_elements()), k(dm3.begin_elements()) ; 
                    i != i_end ; ++i, ++j, ++k)
                {
                    *i = sign_1 * (i.index() % num_limit);
                    *j = sign_2 * (num_limit - (i.index() % num_limit)- 1);
                    *k = sign_1 * sign_2 * ((i.index() % num_limit) * (num_limit - (i.index() % num_limit) - 1));
                    sign_1 *= -1;
                    sign_2 *= (-1) * sign_1;
                }

                ElementProduct<Tag_>::value(dm1, dm2);
                TEST(dm1.lock(lm_read_only),
                        for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()),
                            i_end(dm1.end_elements()),k(dm3.begin_elements()) ; i != i_end ; ++i, ++k)
                        {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, *k, std::numeric_limits<DataType_>::epsilon());
                        },
                        dm1.unlock(lm_read_only));
            }
            DenseMatrix<DataType_> dm01(3, 2, static_cast<DataType_>(1)), dm02(4, 3, static_cast<DataType_>(1)),
                dm03(4, 2, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementProductTest<tags::CPU, float> dense_matrix_elementwise_product_test_float("float");
DenseMatrixElementProductTest<tags::CPU, double> dense_matrix_elementwise_product_test_double("double");
DenseMatrixElementProductTest<tags::CPU::MultiCore, float> mc_dense_matrix_elementwise_product_test_float("MC float");
DenseMatrixElementProductTest<tags::CPU::MultiCore, double> mc_dense_matrix_elementwise_product_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixElementProductTest<tags::CPU::SSE, float> sse_dense_matrix_elementwise_product_test_float("float");
DenseMatrixElementProductTest<tags::CPU::SSE, double> sse_dense_matrix_elementwise_product_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseMatrixElementProductTest<tags::GPU::CUDA, float> cuda_dense_matrix_elementwise_product_test_float("float");
#endif
#ifdef HONEI_CELL
DenseMatrixElementProductTest<tags::Cell, float> cell_dense_matrix_element_product_test_float("Cell float");
DenseMatrixElementProductTest<tags::Cell, double> cell_dense_matrix_element_product_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(0)), dm2(size+1, size, DataType_(0)),
                dm3(size+1, size, DataType_(0));

            unsigned long num_limit(311); //value of used elements will be <= (num_limit/2)^2

            DataType_ sign_1(1), sign_2(1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), 
                i_end(dm1.end_elements()), j(dm2.begin_elements()), k(dm3.begin_elements()) ; 
                i != i_end ; ++i, ++j, ++k)
            {
                *i = sign_1 * (i.index() % num_limit);
                *j = sign_2 * (num_limit - (i.index() % num_limit) - 1);
                *k = sign_1 * sign_2 * (i.index() * (num_limit - (i.index() % num_limit) - 1));
                sign_1 *= -1;
                sign_2 *= (-1) * sign_1;
            }

            ElementProduct<Tag_>::value(dm1, dm2);

            TEST(dm1.lock(lm_read_only), TEST_CHECK_EQUAL(dm1, dm3), dm1.unlock(lm_read_only));

            DenseMatrix<DataType_> dm01(3, 2, static_cast<DataType_>(1)), dm02(4, 3, DataType_(1)),
                dm03(4, 2, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementProductQuickTest<tags::CPU, float> dense_matrix_elementwise_product_quick_test_float("float");
DenseMatrixElementProductQuickTest<tags::CPU, double> dense_matrix_elementwise_product_quick_test_double("double");
DenseMatrixElementProductQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_elementwise_product_quick_test_float("MC float");
DenseMatrixElementProductQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_elementwise_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixElementProductQuickTest<tags::CPU::SSE, float> sse_dense_matrix_elementwise_product_quick_test_float("float");
DenseMatrixElementProductQuickTest<tags::CPU::SSE, double> sse_dense_matrix_elementwise_product_quick_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseMatrixElementProductQuickTest<tags::GPU::CUDA, float> cuda_dense_matrix_elementwise_product_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseMatrixElementProductQuickTest<tags::Cell, float> cell_dense_matrix_element_product_quick_test_float("Cell float");
DenseMatrixElementProductQuickTest<tags::Cell, double> cell_dense_matrix_element_product_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixDenseMatrixElementProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseMatrixElementProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_matrix_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
                SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), 
                    j(sm2.begin_elements()), k(sm3.begin_elements()) ; j != j_end ; ++j, ++k)
                {
                    if (j.index() % 7 == 0)
                    {
                        *j = DataType_(3);
                        *k = DataType_(6);
                    }
                }

                ElementProduct<Tag_>::value(sm2, dm1);

                TEST_CHECK_EQUAL(sm2, sm3);
            }

            DenseMatrix<DataType_> dm03(6, 5);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm01, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm02, dm03), MatrixColumnsDoNotMatch);

            DataType_ factor_1(2.000012345);
            DataType_ factor_2(- 3.000012345);
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(factor_1));
                SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), 
                    j(sm2.begin_elements()), k(sm3.begin_elements()) ; j != j_end ; ++j, ++k)
                {
                    if (j.index() % 7 == 0)
                    {
                        *j = DataType_(factor_2);
                        *k = DataType_(factor_1 * factor_2);
                    }
                }
                ///Test critical elements with operator[] access:
                sm2[0][0] = factor_2;
                sm2[size][size - 1] = factor_2;

                dm1[0][0] = DataType_(factor_1);
                dm1[size][size - 1] = DataType_(factor_1);
                sm3[0][0] = DataType_(factor_1 * factor_2);
                sm3[size][size - 1] = DataType_(factor_1 * factor_2);


                ElementProduct<Tag_>::value(sm2, dm1);

                TEST_CHECK_EQUAL(sm2, sm3);
            }
        }
};
SparseMatrixDenseMatrixElementProductTest<tags::CPU, float> sparse_matrix_dense_matrix_elementwise_product_test_float("float");
SparseMatrixDenseMatrixElementProductTest<tags::CPU, double> sparse_matrix_dense_matrix_elementwise_product_test_double("double");
SparseMatrixDenseMatrixElementProductTest<tags::CPU::MultiCore, float> mc_sparse_matrix_dense_matrix_elementwise_product_test_float("MC float");
SparseMatrixDenseMatrixElementProductTest<tags::CPU::MultiCore, double> mc_sparse_matrix_dense_matrix_elementwise_product_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseMatrixDenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
            SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), j(sm2.begin_elements())                , k(sm3.begin_elements()) ; j != j_end ; ++j, ++k)
            {
                if (j.index() % 7 == 0)
                {
                    *j = DataType_(3);
                    *k = DataType_(6);
                }
            }

            ElementProduct<Tag_>::value(sm2, dm1);

            TEST_CHECK_EQUAL(sm2, sm3);

            DenseMatrix<DataType_> dm03(6, 5);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm01, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm02, dm03), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDenseMatrixElementProductQuickTest<tags::CPU, float> sparse_matrix_dense_matrix_elementwise_product_quick_test_float("float");
SparseMatrixDenseMatrixElementProductQuickTest<tags::CPU, double> sparse_matrix_dense_matrix_elementwise_product_quick_test_double("double");
SparseMatrixDenseMatrixElementProductQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_dense_matrix_elementwise_product_quick_test_float("MC float");
SparseMatrixDenseMatrixElementProductQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_dense_matrix_elementwise_product_quick_test_double("MC double");

template <typename DataType_>
class SparseMatrixBandedMatrixElementProductTest :
    public BaseTest
{
    public:
        SparseMatrixBandedMatrixElementProductTest(const std::string & type) :
            BaseTest("sparse_matrix_banded_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv2(size, DataType_(3));
                BandedMatrix<DataType_> bm2(size, dv2);
                SparseMatrix<DataType_> sm1(size, size, size / 7 + 1), sm3(size, size, size / 7 + 1);
                typename BandedMatrix<DataType_>::ConstElementIterator j(bm2.begin_elements());
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), k(sm3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 7 == 0)
                    {
                        *i = DataType_(2);
                        if (*j != DataType_(0))
                            *k = DataType_(6);
                    }
                }

                ElementProduct<>::value(sm1, bm2);

                TEST_CHECK_EQUAL(sm1, sm3);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, bm01), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixBandedMatrixElementProductTest<float> sparse_matrix_banded_matrix_elementwise_product_test_float("float");
SparseMatrixBandedMatrixElementProductTest<double> sparse_matrix_banded_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class SparseMatrixBandedMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixBandedMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_banded_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> dv2(size, DataType_(3));
            BandedMatrix<DataType_> bm2(size, dv2);
            SparseMatrix<DataType_> sm1(size, size, size / 7 + 1), sm3(size, size, size / 7 + 1);
            typename BandedMatrix<DataType_>::ConstElementIterator j(bm2.begin_elements());
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), k(sm3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 7 == 0)
                {
                    *i = DataType_(2);
                    if (*j != DataType_(0))
                        *k = DataType_(6);
                }
            }

            ElementProduct<>::value(sm1, bm2);

            TEST_CHECK_EQUAL(sm1, sm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, bm01), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixBandedMatrixElementProductQuickTest<float> sparse_matrix_banded_matrix_elementwise_product_quick_test_float("float");
SparseMatrixBandedMatrixElementProductQuickTest<double> sparse_matrix_banded_matrix_elementwise_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class SparseMatrixElementProductTest :
    public BaseTest
{
    public:
        SparseMatrixElementProductTest(const std::string & type) :
            BaseTest("sparse_matrix_elementwise_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                    sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 8 + 1 );
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(2);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_(3);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = DataType_(6);
                    }
                }

                ElementProduct<Tag_>::value(sm1, sm2);

                TEST_CHECK_EQUAL(sm1, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixElementProductTest<tags::CPU, float> sparse_matrix_elementwise_product_test_float("float");
SparseMatrixElementProductTest<tags::CPU, double> sparse_matrix_elementwise_product_test_double("double");
SparseMatrixElementProductTest<tags::CPU::MultiCore, float> mc_sparse_matrix_elementwise_product_test_float("MC float");
SparseMatrixElementProductTest<tags::CPU::MultiCore, double> mc_sparse_matrix_elementwise_product_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(11);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(2);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_(3);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    *k = DataType_(6);
                }
            }

            ElementProduct<Tag_>::value(sm1, sm2);

            TEST_CHECK_EQUAL(sm1, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixElementProductQuickTest<tags::CPU, float> sparse_matrix_elementwise_product_quick_test_float("float");
SparseMatrixElementProductQuickTest<tags::CPU, double> sparse_matrix_elementwise_product_quick_test_double("double");
SparseMatrixElementProductQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_elementwise_product_quick_test_float("MC float");
SparseMatrixElementProductQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_elementwise_product_quick_test_double("MC double");
