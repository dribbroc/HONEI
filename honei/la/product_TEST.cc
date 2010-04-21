/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/product.hh>
#include <honei/la/matrix_error.cc>
#include <honei/la/reduction.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_vector_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0));
                DenseVector<DataType_> dv4(size, DataType_(0));
                DenseVector<DataType_> dv5(dv4.copy());
                for (unsigned long i(0); i < size; ++i)
                    {
                        (dv1)[i]= (-1)*DataType_(i+1);
                        (dv2)[i]= DataType_(size-i);
                        (dv4)[i]= DataType_(i+1);
                        (dv5)[i]= DataType_(i+1);
                    }
                BandedMatrix<DataType_> bm1(size, dv1);
                bm1.insert_band(1, dv4);
                bm1.insert_band(-1, dv5);
                //std::cout << "BM: bm1 = " << bm1 << std::endl;
                //std::cout << "DV: dv = " << dv2 << std::endl;

                DenseVector<DataType_> dv3(size, DataType_(0));
                (dv3)[0]= DataType_(-1);
                (dv3)[size-1]= DataType_(size);
                for (unsigned long i(1); i < size-1; ++i)
                {
                    (dv3)[i]= DataType_((i+1)*(size-i));
                }

                DenseVector<DataType_> prod (Product<Tag_>::value(bm1, dv2));
                //std::cout << "prod = " << prod << std::endl;

                TEST(prod.lock(lm_read_only), TEST_CHECK_EQUAL(prod, dv3), prod.unlock(lm_read_only));
            }

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixDenseVectorProductTest<tags::CPU, float> banded_matrix_dense_vector_product_test_float("float");
BandedMatrixDenseVectorProductTest<tags::CPU, double> banded_matrix_dense_vector_product_test_double("double");
BandedMatrixDenseVectorProductTest<tags::CPU::MultiCore, float> mc_banded_matrix_dense_vector_product_test_float("MC float");
BandedMatrixDenseVectorProductTest<tags::CPU::MultiCore, double> mc_banded_matrix_dense_vector_product_test_double("MC double");
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductTest<tags::CPU::SSE, float> sse_banded_matrix_dense_vector_product_test_float("SSE float");
BandedMatrixDenseVectorProductTest<tags::CPU::SSE, double> sse_banded_matrix_dense_vector_product_test_double("SSE double");
BandedMatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_banded_matrix_dense_vector_product_test_float("MC::SSE float");
BandedMatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_banded_matrix_dense_vector_product_test_double("MC::SSE double");
#endif
#ifdef HONEI_CUDA
BandedMatrixDenseVectorProductTest<tags::GPU::CUDA, float> cuda_banded_matrix_dense_vector_product_test_float("float");
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductTest<tags::Cell, float> cell_banded_matrix_dense_vector_product_test_float("CELL float");
BandedMatrixDenseVectorProductTest<tags::Cell, double> cell_banded_matrix_dense_vector_product_test_double("CELL double");
#endif

template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_vector_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
#if defined HONEI_CELL
            unsigned long size(18225);
#elif defined HONEI_SSE
            unsigned long size(18225);
#else
            unsigned long size(19);
#endif
            unsigned long num_limit(311); //value of used elements will be <= num_limit* size
            DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(1));
            DenseVector<DataType_> dv4(size, DataType_(0));
            DenseVector<DataType_> dv5(dv4.copy());
            DenseVector<DataType_> dv6(size, DataType_(-1)), dv7(size, DataType_(1));
            for (unsigned long i(0); i < size; ++i)
                {
                    (dv1)[i]= (-1)*DataType_((i+1) % num_limit);
                    (dv2)[i]= DataType_((size-i) % num_limit);
                    (dv4)[i]= DataType_((i+1) % num_limit);
                    (dv5)[i]= DataType_((i+1) % num_limit);
                }
            BandedMatrix<DataType_> bm1(size, dv1);
            bm1.insert_band(1, dv4);
            bm1.insert_band(-1, dv5);
            bm1.insert_band(-size + 1, dv6);
            bm1.insert_band(size - 1, dv7);

            DenseVector<DataType_> dv3(size, DataType_(0));
            for (unsigned long i(1); i < size-1; ++i)
            {
                (dv3)[i]= DataType_((dv2[i-1] - dv2[i] + dv2[i+1])*((i+1) % num_limit));
            }
            DenseVector<DataType_> prod(Product<Tag_>::value(bm1, dv2));

            TEST(prod.lock(lm_read_only),
                    for (typename DenseVector<DataType_>::ConstElementIterator dit(dv3.begin_elements()), it(prod.begin_elements()), i_end(prod.end_elements()) ;
                        it != i_end ; ++it, ++dit)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*it, *dit, std::numeric_limits<DataType_>::epsilon());
                    },
                    prod.unlock(lm_read_only));

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixDenseVectorProductQuickTest<tags::CPU, float> banded_matrix_dense_vector_product_quick_test_float("float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU, double> banded_matrix_dense_vector_product_quick_test_double("double");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, float> mc_banded_matrix_dense_vector_product_quick_test_float("MC float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, double> mc_banded_matrix_dense_vector_product_quick_test_double("MC double");
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductQuickTest<tags::CPU::SSE, float> sse_banded_matrix_dense_vector_product_quick_test_float("SSE float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_banded_matrix_dense_vector_product_quick_test_double("SSE double");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_banded_matrix_dense_vector_product_quick_test_float("MC::SSE float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_banded_matrix_dense_vector_product_quick_test_double("MC::SSE double");
#endif
#ifdef HONEI_CUDA
BandedMatrixDenseVectorProductQuickTest<tags::GPU::CUDA, float> cuda_banded_matrix_dense_vector_product_quick_test_float("float");
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductQuickTest<tags::Cell, float> cell_banded_matrix_dense_vector_product_quick_test_float("CELL float");
BandedMatrixDenseVectorProductQuickTest<tags::Cell, double> cell_banded_matrix_dense_vector_product_quick_test_double("CELL double");
#endif

template <typename Tag_, typename DataType_>
class Q1MatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        Q1MatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("q1_matrix_dense_vector_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size;
            for (unsigned long level(1) ; level <= 10 ; ++level)
            {
                size  = (unsigned long)pow((pow(2, level) + 1), 2);
                unsigned long num_limit(311); //value of used elements will be <= num_limit* size

                DenseVector<DataType_> dv1(size, DataType_(1));
                DenseVector<DataType_> dv2(size, DataType_(1));
                DenseVector<DataType_> dv3(size, DataType_(1));
                DenseVector<DataType_> dv4(size, DataType_(1));

                for (unsigned long i(0); i < size; ++i)
                {
                    (dv1)[i]= DataType_((i + 2) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv2)[i]= DataType_((i + 1) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv3)[i]= DataType_((i + 25) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv4)[i]= DataType_((i + 7) % num_limit);
                }

                BandedMatrixQ1<DataType_> bm1(size, dv3, dv2, dv4, dv2, dv3, dv4, dv4, dv3, dv2);
                BandedMatrix<DataType_> bm2(bm1);

                DenseVector<DataType_> prod(Product<Tag_>::value(bm1, dv1));
#ifdef HONEI_SSE
            DenseVector<DataType_> dv_ref(Product<tags::CPU::SSE>::value(bm2, dv1));
#else
            DenseVector<DataType_> dv_ref(Product<tags::CPU>::value(bm2, dv1));
#endif

            TEST(prod.lock(lm_read_only),
                    for (typename DenseVector<DataType_>::ConstElementIterator dit(dv_ref.begin_elements()), it(prod.begin_elements()), i_end(prod.end_elements()) ;
                        it != i_end ; ++it, ++dit)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*it, *dit, std::numeric_limits<DataType_>::epsilon());
                    },
                    prod.unlock(lm_read_only));
            }

            DenseVector<DataType_> dv01(4, DataType_(1));
            DenseVector<DataType_> dv02(1089, DataType_(1));
            BandedMatrixQ1<DataType_> bm01(1089, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02);
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
Q1MatrixDenseVectorProductTest<tags::CPU, float> q1_prod_test_float("float");
Q1MatrixDenseVectorProductTest<tags::CPU, double> q1_prod_test_double("double");
Q1MatrixDenseVectorProductTest<tags::CPU::MultiCore, float> q1_prod_mc_test_float("MC float");
Q1MatrixDenseVectorProductTest<tags::CPU::MultiCore, double> q1_prod_mc_test_double("MC double");
#ifdef HONEI_SSE
Q1MatrixDenseVectorProductTest<tags::CPU::SSE, float> sse_q1_prod_test_float("float");
Q1MatrixDenseVectorProductTest<tags::CPU::SSE, double> sse_q1_prod_test_double("double");
Q1MatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, float> q1_prod_mc_sse_test_float("MC SSE float");
Q1MatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, double> q1_prod_mc_sse_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
Q1MatrixDenseVectorProductTest<tags::GPU::CUDA, float> cuda_q1_prod_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixDenseVectorProductTest<tags::GPU::CUDA, double> cuda_q1_prod_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class Q1MatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        Q1MatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("q1_matrix_dense_vector_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(1025ul * 1025);
            unsigned long num_limit(311); //value of used elements will be <= num_limit* size

            DenseVector<DataType_> dv1(size, DataType_(1));
            DenseVector<DataType_> dv2(size, DataType_(1));
            DenseVector<DataType_> dv3(size, DataType_(1));
            DenseVector<DataType_> dv4(size, DataType_(1));

            for (unsigned long i(0); i < size; ++i)
            {
                (dv1)[i]= DataType_((i + 2) % num_limit);
            }
            for (unsigned long i(0) ; i < size ; ++i)
            {
                (dv2)[i]= DataType_((i + 1) % num_limit);
            }
            for (unsigned long i(0) ; i < size ; ++i)
            {
                (dv3)[i]= DataType_((i + 25) % num_limit);
            }
            for (unsigned long i(0) ; i < size ; ++i)
            {
                (dv4)[i]= DataType_((i + 7) % num_limit);
            }

            BandedMatrixQ1<DataType_> bm1(size, dv3, dv2, dv4, dv2, dv3, dv4, dv4, dv3, dv2);
            BandedMatrix<DataType_> bm2(bm1);

            DenseVector<DataType_> prod(Product<Tag_>::value(bm1, dv1));
#ifdef HONEI_SSE
            DenseVector<DataType_> dv_ref(Product<tags::CPU::SSE>::value(bm2, dv1));
#else
            DenseVector<DataType_> dv_ref(Product<tags::CPU>::value(bm2, dv1));
#endif

            TEST(prod.lock(lm_read_only),
                    for (typename DenseVector<DataType_>::ConstElementIterator dit(dv_ref.begin_elements()), it(prod.begin_elements()), i_end(prod.end_elements()) ;
                        it != i_end ; ++it, ++dit)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*it, *dit, std::numeric_limits<DataType_>::epsilon());
                    },
                    prod.unlock(lm_read_only));

            DenseVector<DataType_> dv01(4, DataType_(1));
            DenseVector<DataType_> dv02(1089, DataType_(1));
            BandedMatrixQ1<DataType_> bm01(1089, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02);
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
Q1MatrixDenseVectorProductQuickTest<tags::CPU, float> q1_prod_quick_test_float("float");
Q1MatrixDenseVectorProductQuickTest<tags::CPU, double> q1_prod_quick_test_double("double");
Q1MatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, float> q1_prod_quick_mc_test_float("MC float");
Q1MatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, double> q1_prod_quick_mc_test_double("MC double");
#ifdef HONEI_SSE
Q1MatrixDenseVectorProductQuickTest<tags::CPU::SSE, float> sse_q1_prod_quick_test_float("float");
Q1MatrixDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_q1_prod_quick_test_double("double");
Q1MatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, float> q1_prod_quick_mc_sse_test_float("MC SSE float");
Q1MatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, double> q1_prod_quick_mc_sse_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
Q1MatrixDenseVectorProductQuickTest<tags::GPU::CUDA, float> cuda_q1_prod_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixDenseVectorProductQuickTest<tags::GPU::CUDA, double> cuda_q1_prod_quick_test_double("double");
#endif
#endif

template <typename DataType_>
class BandedMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        BandedMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseVector<DataType_> dv4(size, DataType_(2));
                DenseVector<DataType_> dv5(dv4.copy());
                bm1.insert_band(1, dv4);
                bm1.insert_band(-1, dv5);
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(6 * 1);
                    if (i.index() % 10 == 1) *i = DataType_(6 * 1);
                    if (i.index() % 10 == 10 - 1) *i = DataType_(6 * 1);
                }
                if (!(size % 12 == 0))
                {
                   sv2[size-1] = DataType_(0);
                }

                SparseVector<DataType_> prod(Product<>::value(bm1, sv1));

                TEST_CHECK_EQUAL(prod, sv2);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseVector<DataType_> sv01(4, 3);
            TEST_CHECK_THROWS(Product<>::value(bm01, sv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixSparseVectorProductTest<float> banded_matrix_sparse_vector_product_test_float("float");
BandedMatrixSparseVectorProductTest<double> banded_matrix_sparse_vector_product_test_double("double");

template <typename DataType_>
class BandedMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> dv1(size, DataType_(2));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseVector<DataType_> dv4(size, DataType_(2));
            DenseVector<DataType_> dv5(dv4.copy());
            bm1.insert_band(1, dv4);
            bm1.insert_band(-1, dv5);
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 5 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 5 == 0) *i = DataType_(6 * 1);
                if (i.index() % 5 == 1) *i = DataType_(6 * 1);
                if (i.index() % 5 == 5 - 1) *i = DataType_(6 * 1);
            }

            sv2[size-1] = DataType_(0);
            SparseVector<DataType_> prod(Product<>::value(bm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

            BandedMatrix<DataType_> bm01(5);
            SparseVector<DataType_> sv01(4, 3);
            TEST_CHECK_THROWS(Product<>::value(bm01, sv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixSparseVectorProductQuickTest<float> banded_matrix_sparse_vector_product_quick_test_float("float");
BandedMatrixSparseVectorProductQuickTest<double> banded_matrix_sparse_vector_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class DenseMatrixDenseVectorRangeProductTest :
    public BaseTest
{
    public:
        DenseMatrixDenseVectorRangeProductTest(const std::string & type) :
            BaseTest("dense_matrix_dense_vector_range_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 11) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
                DenseVector<DataType_> dv1(2 * size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
                DenseVectorRange<DataType_> dvr1(dv1, size, 0);
                DenseVector<DataType_> prod(Product<Tag_>::value(dm1, dvr1));

                TEST_CHECK_EQUAL(prod, dv2);
            }

            DenseMatrix<DataType_> dm01(4, 3, DataType_(1));
            DenseVector<DataType_> dv01(4, DataType_(1));
            DenseVectorRange<DataType_> dvr01(dv01, 4, 0);

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dvr01), VectorSizeDoesNotMatch);
        }
};

DenseMatrixDenseVectorRangeProductTest<tags::CPU, float> dense_matrix_dense_vector_r_product_test_float("float");
DenseMatrixDenseVectorRangeProductTest<tags::CPU, double> dense_matrix_dense_vector_r_product_test_double("double");
DenseMatrixDenseVectorRangeProductTest<tags::CPU::MultiCore, float> mc_dense_matrix_dense_vector_r_product_test_float("MC float");
DenseMatrixDenseVectorRangeProductTest<tags::CPU::MultiCore, double> mc_dense_matrix_dense_vector_r_product_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixDenseVectorRangeProductTest<tags::CPU::SSE, float> sse_dense_matrix_dense_vector_product_r_test_float("SSE float");
DenseMatrixDenseVectorRangeProductTest<tags::CPU::SSE, double> sse_dense_matrix_dense_vector_product_r_test_double("SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixDenseVectorRangeProductTest<tags::Cell, float> cell_dense_matrix_dense_vector_product_r_test_float("Cell float");
DenseMatrixDenseVectorRangeProductTest<tags::Cell, double> cell_dense_matrix_dense_vector_product_r_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        DenseMatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("dense_matrix_dense_vector_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 11) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
                DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
                DenseVector<DataType_> prod(Product<Tag_>::value(dm1, dv1));

                TEST_CHECK_EQUAL(prod, dv2);
            }

            DenseMatrix<DataType_> dm01(4, 3, DataType_(1));
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dv01), VectorSizeDoesNotMatch);
        }
};

DenseMatrixDenseVectorProductTest<tags::CPU, float> dense_matrix_dense_vector_product_test_float("float");
DenseMatrixDenseVectorProductTest<tags::CPU, double> dense_matrix_dense_vector_product_test_double("double");
// MultiCore
DenseMatrixDenseVectorProductTest<tags::CPU::MultiCore, float> mc_dense_matrix_dense_vector_product_test_float("MC float");
DenseMatrixDenseVectorProductTest<tags::CPU::MultiCore, double> mc_dense_matrix_dense_vector_product_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixDenseVectorProductTest<tags::CPU::SSE, float> sse_dense_matrix_dense_vector_product_test_float("SSE float");
DenseMatrixDenseVectorProductTest<tags::CPU::SSE, double> sse_dense_matrix_dense_vector_product_test_double("SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixDenseVectorProductTest<tags::Cell, float> cell_dense_matrix_dense_vector_product_test_float("Cell float");
DenseMatrixDenseVectorProductTest<tags::Cell, double> cell_dense_matrix_dense_vector_product_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_dense_vector_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(31);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
            DenseVector<DataType_> prod(Product<Tag_>::value(dm1, dv1));

            TEST_CHECK_EQUAL(prod, dv2);

            DenseMatrix<DataType_> dm01(4, 3, DataType_(1));
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dv01), VectorSizeDoesNotMatch);

        }
};
DenseMatrixDenseVectorProductQuickTest<tags::CPU, float> dense_matrix_dense_vector_product_quick_test_float("float");
DenseMatrixDenseVectorProductQuickTest<tags::CPU, double> dense_matrix_dense_vector_product_quick_test_double("double");
// MultiCore
DenseMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_dense_vector_product_quick_test_float("MC float");
DenseMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_dense_vector_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixDenseVectorProductQuickTest<tags::CPU::SSE, float> sse_dense_matrix_dense_vector_product_quick_test_float("SSE float");
DenseMatrixDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_dense_matrix_dense_vector_product_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixDenseVectorProductQuickTest<tags::Cell, float> cell_dense_matrix_dense_vector_product_quick_test_float("Cell float");
DenseMatrixDenseVectorProductQuickTest<tags::Cell, double> cell_dense_matrix_dense_vector_product_quick_test_double("Cell double");
#endif

template <typename DataType_>
class DenseMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        DenseMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(6 * (size / 10 + 1));
                }
                DenseVector<DataType_> prod(Product<>::value(dm1, sv1));

                for (typename DenseVector<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()) ;
                    i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(6 * (size / 10 + 1)), std::numeric_limits<DataType_>::epsilon());
                }
            }

            DenseMatrix<DataType_> dm01(4, 3, DataType_(1));
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(Product<>::value(dm01, sv01), VectorSizeDoesNotMatch);
        }
};

DenseMatrixSparseVectorProductTest<float> dense_matrix_sparse_vector_product_test_float("float");
DenseMatrixSparseVectorProductTest<double> dense_matrix_sparse_vector_product_test_double("double");

template <typename DataType_>
class DenseMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            DenseVector<DataType_> prod(Product<>::value(dm1, sv1));
            for (typename DenseVector<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()) ;
                i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(12), std::numeric_limits<DataType_>::epsilon());
            }

            DenseMatrix<DataType_> dm01(4, 3, DataType_(1));
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(Product<>::value(dm01, sv01), VectorSizeDoesNotMatch);
        }
};

DenseMatrixSparseVectorProductQuickTest<float> dense_matrix_sparse_vector_product_test_quick_float("float");
DenseMatrixSparseVectorProductQuickTest<double> dense_matrix_sparse_vector_product_test_quick_double("double");

template <typename DataType_>
class SparseMatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size,  size / 8 + 1);
                DenseVector<DataType_> dv1(size,  DataType_(3)), dv2(size + 1, DataType_(6 * size));
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                DenseVector<DataType_> prod(Product<>::value(sm1, dv1));

                for (typename DenseVector<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()) ;
                    i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(6 * size), std::numeric_limits<DataType_>::epsilon());
                }
            }

            SparseMatrix<DataType_> sm01(4, 3, 1);
            DenseVector<DataType_> dv01(4, DataType_(3));

            TEST_CHECK_THROWS(Product<>::value(sm01, dv01), VectorSizeDoesNotMatch);
        }
};
SparseMatrixDenseVectorProductTest<float> sparse_matrix_dense_vector_product_test_float("float");
SparseMatrixDenseVectorProductTest<double> sparse_matrix_dense_vector_product_test_double("double");

template <typename DataType_>
class SparseMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
            DenseVector<DataType_> prod(Product<>::value(sm1, dv1));

            for (typename DenseVector<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()) ;
                    i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(6 * size), std::numeric_limits<DataType_>::epsilon());
            }

            SparseMatrix<DataType_> sm01(4, 3, 1);
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(Product<>::value(sm01, dv01), VectorSizeDoesNotMatch);

        }
};
SparseMatrixDenseVectorProductQuickTest<float> sparse_matrix_dense_vector_product_quick_test_float("float");
SparseMatrixDenseVectorProductQuickTest<double> sparse_matrix_dense_vector_product_quick_test_double("double");

template <typename DataType_, typename Tag_>
class SparseMatrixELLDenseVectorProductTest :
    public BaseTest
{
    public:
        SparseMatrixELLDenseVectorProductTest(const std::string & type) :
            BaseTest("sparse_matrix_ell_dense_vector_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sms(size, size + 3);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sms.begin_elements()) ; i < sms.end_elements() ; ++i)
                {
                    if (i.index() % 5 == 0)
                        *i = DataType_(i.index()) / 1.234;
                }
                DenseMatrix<DataType_> dm0(sms);
                SparseMatrixELL<DataType_> sm0(sms);
                DenseVector<DataType_> dv1(size + 3, DataType_(4));
                DenseVector<DataType_> result(size, DataType_(4711));
                dv1[0] = 1;
                dv1[1] = 2;
                DenseVector<DataType_> prod(Product<Tag_>::value(result, sm0, dv1));
                DenseVector<DataType_> prod_ref(Product<>::value(dm0, dv1));

                prod.lock(lm_read_only);
                TEST_CHECK_EQUAL(prod, prod_ref);
                prod.unlock(lm_read_only);
            }
        }
};
SparseMatrixELLDenseVectorProductTest<float, tags::CPU> sparse_matrix_ell_dense_vector_product_test_float("float");
SparseMatrixELLDenseVectorProductTest<double, tags::CPU> sparse_matrix_ell_dense_vector_product_test_double("double");
#ifdef HONEI_SSE
SparseMatrixELLDenseVectorProductTest<float, tags::CPU::SSE> sse_sparse_matrix_ell_dense_vector_product_test_float("float");
SparseMatrixELLDenseVectorProductTest<double, tags::CPU::SSE> sse_sparse_matrix_ell_dense_vector_product_test_double("double");
SparseMatrixELLDenseVectorProductTest<float, tags::CPU::MultiCore::SSE> mc_sse_sparse_matrix_ell_dense_vector_product_test_float("float");
SparseMatrixELLDenseVectorProductTest<double, tags::CPU::MultiCore::SSE> mc_sse_sparse_matrix_ell_dense_vector_product_test_double("double");
#endif
#ifdef HONEI_OPENCL
SparseMatrixELLDenseVectorProductTest<float, tags::OpenCL::CPU> ocl_cpu_sparse_matrix_ell_dense_vector_product_test_float("float");
SparseMatrixELLDenseVectorProductTest<double, tags::OpenCL::CPU> ocl_cpu_sparse_matrix_ell_dense_vector_product_test_double("double");
SparseMatrixELLDenseVectorProductTest<float, tags::OpenCL::GPU> ocl_gpu_sparse_matrix_ell_dense_vector_product_test_float("float");
SparseMatrixELLDenseVectorProductTest<double, tags::OpenCL::GPU> ocl_gpu_sparse_matrix_ell_dense_vector_product_test_double("double");
#endif
#ifdef HONEI_CUDA
SparseMatrixELLDenseVectorProductTest<float, tags::GPU::CUDA> cuda_sparse_matrix_ell_dense_vector_product_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
SparseMatrixELLDenseVectorProductTest<double, tags::GPU::CUDA> cuda_sparse_matrix_ell_dense_vector_product_test_double("double");
#endif
#endif

template <typename DataType_, typename Tag_>
class SparseMatrixELLDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixELLDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_ell_dense_vector_product_quick_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size (1000);
            SparseMatrix<DataType_> sms(size, size + 3);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sms.begin_elements()) ; i < sms.end_elements() ; ++i)
            {
                if (i.index() % 5 == 0)
                    *i = DataType_(i.index()) / 1.234;
            }
            DenseMatrix<DataType_> dm0(sms);
            SparseMatrixELL<DataType_> sm0(sms);
            DenseVector<DataType_> dv1(size + 3, DataType_(4));
            DenseVector<DataType_> result(size, DataType_(4711));
            dv1[0] = 1;
            dv1[1] = 2;
            DenseVector<DataType_> prod(Product<Tag_>::value(result, sm0, dv1));
            DenseVector<DataType_> prod_ref(Product<tags::CPU>::value(dm0, dv1));

            prod.lock(lm_read_only);
            TEST_CHECK_EQUAL(prod, prod_ref);
            prod.unlock(lm_read_only);
        }
};
SparseMatrixELLDenseVectorProductQuickTest<float, tags::CPU> sparse_matrix_ell_dense_vector_product_quick_test_float("float");
SparseMatrixELLDenseVectorProductQuickTest<double, tags::CPU> sparse_matrix_ell_dense_vector_product_quick_test_double("double");
#ifdef HONEI_SSE
SparseMatrixELLDenseVectorProductQuickTest<float, tags::CPU::SSE> sse_sparse_matrix_ell_dense_vector_product_quick_test_float("float");
SparseMatrixELLDenseVectorProductQuickTest<double, tags::CPU::SSE> sse_sparse_matrix_ell_dense_vector_product_quick_test_double("double");
SparseMatrixELLDenseVectorProductQuickTest<float, tags::CPU::MultiCore::SSE> mc_sse_sparse_matrix_ell_dense_vector_product_quick_test_float("float");
SparseMatrixELLDenseVectorProductQuickTest<double, tags::CPU::MultiCore::SSE> mc_sse_sparse_matrix_ell_dense_vector_product_quick_test_double("double");
#endif
#ifdef HONEI_OPENCL
SparseMatrixELLDenseVectorProductQuickTest<float, tags::OpenCL::CPU> ocl_cpu_sparse_matrix_ell_dense_vector_product_quick_test_float("float");
SparseMatrixELLDenseVectorProductQuickTest<double, tags::OpenCL::CPU> ocl_cpu_sparse_matrix_ell_dense_vector_product_quick_test_double("double");
SparseMatrixELLDenseVectorProductQuickTest<float, tags::OpenCL::GPU> ocl_gpu_sparse_matrix_ell_dense_vector_product_quick_test_float("float");
SparseMatrixELLDenseVectorProductQuickTest<double, tags::OpenCL::GPU> ocl_gpu_sparse_matrix_ell_dense_vector_product_quick_test_double("double");
#endif
#ifdef HONEI_CUDA
SparseMatrixELLDenseVectorProductQuickTest<float, tags::GPU::CUDA> cuda_sparse_matrix_ell_dense_vector_product_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
SparseMatrixELLDenseVectorProductQuickTest<double, tags::GPU::CUDA> cuda_sparse_matrix_ell_dense_vector_product_quick_test_double("double");
#endif
#endif

template <typename DataType_>
class SparseMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        SparseMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("sparse_matrix_sparse_vector_product_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                        i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(6 * size);
                }
                SparseVector<DataType_> prod(Product<>::value(sm1, sv1));

                TEST_CHECK_EQUAL(prod, sv2);

                SparseMatrix<DataType_> sm01(4, 3, 1);
                SparseVector<DataType_> sv01(4, 3);

                TEST_CHECK_THROWS(Product<>::value(sm01, sv01), VectorSizeDoesNotMatch);
            }
        }
};
SparseMatrixSparseVectorProductTest<float> sparse_matrix_sparse_vector_product_test_float("float");
SparseMatrixSparseVectorProductTest<double> sparse_matrix_sparse_vector_product_test_double("double");

template <typename DataType_>
class SparseMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_sparse_vector_product_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(20);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(12);
            }
            SparseVector<DataType_> prod(Product<>::value(sm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

            SparseMatrix<DataType_> sm01(4, 3, 1);
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(Product<>::value(sm01, sv01), VectorSizeDoesNotMatch);
        }
};
SparseMatrixSparseVectorProductQuickTest<float> sparse_matrix_sparse_vector_product_test_quick_float("float");
SparseMatrixSparseVectorProductQuickTest<double> sparse_matrix_sparse_vector_product_test_quick_double("double");

template <typename DataType_>
class BandedMatrixProductTest :
    public BaseTest
{
    public:
        BandedMatrixProductTest(const std::string & type) :
            BaseTest("banded_matrix_product_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                DenseVector<DataType_> dv3(size, DataType_(6));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);
                BandedMatrix<DataType_> prod(Product<>::value(bm1, bm2));

                TEST_CHECK_EQUAL(prod, bm3);
            }
            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Product<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixProductTest<float> banded_matrix_product_test_float("float");
BandedMatrixProductTest<double> banded_matrix_product_test_double("double");

template <typename DataType_>
class BandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_product_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            DenseVector<DataType_> dv4(size, DataType_(6));

            BandedMatrix<DataType_> bm1(size, dv1);
            bm1.insert_band(1, dv3.copy());
            bm1.insert_band(-2, dv4.copy());

            DenseVector<DataType_> dv2(size, DataType_(3));
            BandedMatrix<DataType_> bm2(size, dv2.copy());
            bm2.insert_band(size - 1, dv2.copy());
            bm2.insert_band(1 - size, dv2.copy());

            //SparseMatrix equal to bm2
            SparseMatrix<DataType_> sm1(size, size);
            sm1[0][size-1] = DataType_(3);
            sm1[size-1][0] = DataType_(3);
            for(unsigned long i=0; i < size ; ++i)
            {
                sm1(i, i) = DataType_(3);
            }

            //Calculate products for banded * banded and banded * sparse and check equality.
            BandedMatrix<DataType_> prod(Product<>::value(bm1, bm2));
            DenseMatrix<DataType_> prod2(Product<>::value(bm1, sm1));
            typename BandedMatrix<DataType_>::ConstElementIterator k(prod.begin_elements());
            for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod2.begin_elements()), l_end(prod2.end_elements()) ; l != l_end ; ++l, ++k)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Product<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixProductQuickTest<float> banded_matrix_product_quick_test_float("float");
BandedMatrixProductQuickTest<double> banded_matrix_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class DenseMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_product_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long i(0);
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1, ++i)
            {
                DenseMatrix<DataType_> dm1(size + 3, size + 2);
                DenseMatrix<DataType_> dm2(size + 2, size + 1);
                DenseMatrix<DataType_> dm3(size + 3 , size + 1, DataType_(0));
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 1.23456;
                }
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 3.9876;
                }

                typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
                for (unsigned int s(0) ; s < dm1.rows() ; ++s)
                {
                    const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                    for (unsigned int t(0); t < dm2.columns() ; ++t)
                    {
                        const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                        *i = DotProduct<>::value(dm2_column, dm1_row);
                        ++i;
                    }
                }

                DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));
                for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                        j(dm3.begin_elements()); i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i / *j, 1, std::numeric_limits<DataType_>::epsilon() * 3);
                }
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1)), dm02(3, 3, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dm02), MatrixRowsDoNotMatch);
        }
};
DenseMatrixProductTest<tags::CPU, float> dense_matrix_product_test_float("float");
DenseMatrixProductTest<tags::CPU, double> dense_matrix_product_test_double("double");
DenseMatrixProductTest<tags::CPU::MultiCore, float> mc_dense_matrix_product_test_float("MC float");
DenseMatrixProductTest<tags::CPU::MultiCore, double> mc_dense_matrix_product_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixProductTest<tags::CPU::SSE, float> sse_dense_matrix_product_test_float("SSE float");
DenseMatrixProductTest<tags::CPU::SSE, double> sse_dense_matrix_product_test_double("SSE double");
DenseMatrixProductTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_product_test_float("MC SSE float");
DenseMatrixProductTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_product_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_product_quick_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size(32);

            DenseMatrix<DataType_> dm1(size + 3, size + 2);
            DenseMatrix<DataType_> dm2(size + 2, size + 1);
            DenseMatrix<DataType_> dm3(size + 3 , size + 1, DataType_(0));
            for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(i.index()) / 1.23456;
            }
            for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(i.index()) / 3.9876;
            }

            typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
            for (unsigned int s(0) ; s < dm1.rows() ; ++s)
            {
                const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                for (unsigned int t(0); t < dm2.columns() ; ++t)
                {
                    const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                    *i = DotProduct<>::value(dm2_column, dm1_row);
                    ++i;
                }
            }

            DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));
            for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                    j(dm3.begin_elements()); i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i / *j, 1, std::numeric_limits<DataType_>::epsilon() * 3);
            }

            DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(3, 4, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm02, dm01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixProductQuickTest<tags::CPU, float>  dense_matrix_product_quick_test_float("float");
DenseMatrixProductQuickTest<tags::CPU, double> dense_matrix_product_quick_test_double("double");
DenseMatrixProductQuickTest<tags::CPU::MultiCore, float>  mc_dense_matrix_product_quick_test_float("MC float");
DenseMatrixProductQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixProductQuickTest<tags::CPU::SSE, float>  sse_dense_matrix_product_quick_test_float("SSE float");
DenseMatrixProductQuickTest<tags::CPU::SSE, double> sse_dense_matrix_product_quick_test_double("SSE double");
DenseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, float>  sse_mc_dense_matrix_product_quick_test_float("MC SSE float");
DenseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_product_quick_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductNX2Test :
    public BaseTest
{
    public:
        DenseMatrixProductNX2Test(const std::string & type) :
            BaseTest("dense_matrix_nx2_product_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long i(0);
            for (unsigned long size(10) ; size < (1 << 11) ; size <<= 1, ++i)
            {
                DenseMatrix<DataType_> dm1(size, size);
                DenseMatrix<DataType_> dm2(size, 2);
                DenseMatrix<DataType_> dm3(size , 2, DataType_(0));
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 1.23456;
                }
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 3.9876;
                }

                typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
                for (unsigned int s(0) ; s < dm1.rows() ; ++s)
                {
                    const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                    for (unsigned int t(0); t < dm2.columns() ; ++t)
                    {
                        const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                        *i = DotProduct<>::value(dm2_column, dm1_row);
                        ++i;
                    }
                }

                DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));
                for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                        j(dm3.begin_elements()); i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1)), dm02(3, 3, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dm02), MatrixRowsDoNotMatch);
        }
};
DenseMatrixProductNX2Test<tags::CPU, float> dense_matrix_product_nx2_test_float("float");
DenseMatrixProductNX2Test<tags::CPU, double> dense_matrix_product_nx2_test_double("double");
DenseMatrixProductNX2Test<tags::CPU::MultiCore, float> mc_dense_matrix_product_nx2_test_float("MC float");
DenseMatrixProductNX2Test<tags::CPU::MultiCore, double> mc_dense_matrix_product_nx2_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixProductNX2Test<tags::CPU::SSE, float> sse_dense_matrix_product_nx2_test_float("SSE float");
DenseMatrixProductNX2Test<tags::CPU::SSE, double> sse_dense_matrix_product_nx2_test_double("SSE double");
DenseMatrixProductNX2Test<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_product_nx2_test_float("MC SSE float");
DenseMatrixProductNX2Test<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_product_nx2_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductNX2QuickTest :
    public QuickTest
{
    public:
        DenseMatrixProductNX2QuickTest(const std::string & type) :
            QuickTest("dense_matrix_product_nx2_quick_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size(32);

            DenseMatrix<DataType_> dm1(size, size);
            DenseMatrix<DataType_> dm2(size, 2);
            DenseMatrix<DataType_> dm3(size , 2, DataType_(0));
            for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(i.index()) / 1.23456;
            }
            for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(i.index()) / 3.9876;
            }

            typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
            for (unsigned int s(0) ; s < dm1.rows() ; ++s)
            {
                const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                for (unsigned int t(0); t < dm2.columns() ; ++t)
                {
                    const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                    *i = DotProduct<>::value(dm2_column, dm1_row);
                    ++i;
                }
            }

            DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));
            for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                    j(dm3.begin_elements()); i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
            }

            DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(3, 4, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm02, dm01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixProductNX2QuickTest<tags::CPU, float>  dense_matrix_product_nx2_quick_test_float("float");
DenseMatrixProductNX2QuickTest<tags::CPU, double> dense_matrix_product_nx2_quick_test_double("double");
DenseMatrixProductNX2QuickTest<tags::CPU::MultiCore, float>  mc_dense_matrix_product_nx2_quick_test_float("MC float");
DenseMatrixProductNX2QuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_product_nx2_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixProductNX2QuickTest<tags::CPU::SSE, float>  sse_dense_matrix_product_nx2_quick_test_float("SSE float");
DenseMatrixProductNX2QuickTest<tags::CPU::SSE, double> sse_dense_matrix_product_nx2_quick_test_double("SSE double");
DenseMatrixProductNX2QuickTest<tags::CPU::MultiCore::SSE, float>  sse_mc_dense_matrix_product_nx2_quick_test_float("MC SSE float");
DenseMatrixProductNX2QuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_product_nx2_quick_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductCellTest :
    public BaseTest
{
    public:
        DenseMatrixProductCellTest(const std::string & type) :
            BaseTest("dense_matrix_product_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                for (unsigned i(0) ; i < 4 ; i++)
                {
                    DenseMatrix<DataType_> dm1(size + i + 1, size + 2);
                    DenseMatrix<DataType_> dm2(size + 2, size + i);
                    DenseMatrix<DataType_> dm3(size + i + 1 , size + i, DataType_(0));
                    for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                            i != i_end ; ++i)
                    {
                        *i = DataType_(i.index()) / 1.23456;
                    }
                    for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                            i != i_end ; ++i)
                    {
                        *i = DataType_(i.index()) / 3.9876;
                    }

                    typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
                    for (unsigned int s(0) ; s < dm1.rows() ; ++s)
                    {
                        const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                        for (unsigned int t(0); t < dm2.columns() ; ++t)
                        {
                            const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                            *i = DotProduct<>::value(dm2_column, dm1_row);
                            ++i;
                        }
                    }
                    DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));

                    for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                            j(dm3.begin_elements()); i != i_end ; ++i, ++j)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, (size / 3.5 * *i) * std::numeric_limits<DataType_>::epsilon());
                    }
                }
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1)), dm02(3, 3, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm01, dm02), MatrixRowsDoNotMatch);
        }
};

#ifdef HONEI_CELL
DenseMatrixProductCellTest<tags::Cell, float> cell_dense_matrix_product_test_float("Cell float");
DenseMatrixProductCellTest<tags::Cell, double> cell_dense_matrix_product_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class DenseMatrixProductCellQuickTest :
    public QuickTest
{
    public:
        DenseMatrixProductCellQuickTest(const std::string & type) :
            QuickTest("dense_matrix_product_quick_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size(192);

            for (unsigned i(0) ; i < 4 ; i++)
            {
                DenseMatrix<DataType_> dm1(size + 3, size + 2);
                DenseMatrix<DataType_> dm2(size + 2, size + i);
                DenseMatrix<DataType_> dm3(size + 3 , size + i, DataType_(0));

                for (typename DenseMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 1.23456;
                }
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 3.9876;
                }

                typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
                for (unsigned int s(0) ; s < dm1.rows() ; ++s)
                {
                    const DenseVectorRange<DataType_> dm1_row(dm1[s]);
                    for (unsigned int t(0); t < dm2.columns() ; ++t)
                    {
                        const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                        *i = DotProduct<>::value(dm2_column, dm1_row);
                        ++i;
                    }
                }
                DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));
                for(typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                        j(dm3.begin_elements()); i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, (size / 3.5 * *i) * std::numeric_limits<DataType_>::epsilon());
                }

                DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(3, 4, DataType_(1));

                TEST_CHECK_THROWS(Product<Tag_>::value(dm02, dm01), MatrixRowsDoNotMatch);
            }
        }
};
#ifdef HONEI_CELL
DenseMatrixProductCellQuickTest<tags::Cell, float> cell_dense_matrix_product_quick_test_float("Cell float");
DenseMatrixProductCellQuickTest<tags::Cell, double> cell_dense_matrix_product_quick_test_double("Cell double");
#endif

template <typename DataType_>
class DenseMatrixSparseMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_product_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
                SparseMatrix<DataType_> sm2(size, size+1, size / 8 + 1);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                        i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }
                DenseMatrix<DataType_> prod(Product<>::value(dm1, sm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            SparseMatrix<DataType_> sm01(3, 3, 1);
            DenseMatrix<DataType_> dm02(3, 4);

            TEST_CHECK_THROWS(Product<>::value(dm02, sm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixSparseMatrixProductTest<float> dense_matrix_sparse_matrix_product_test_float("float");
DenseMatrixSparseMatrixProductTest<double> dense_matrix_sparse_matrix_product_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_product_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
            SparseMatrix<DataType_> sm2(size, size+1, size / 8 + 1);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                    i_end(sm2.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }
            DenseMatrix<DataType_> prod(Product<>::value(dm1, sm2));

            TEST_CHECK_EQUAL(prod, dm3);

            SparseMatrix<DataType_> sm01(3, 3, 1);
            DenseMatrix<DataType_> dm02(3, 4);

            TEST_CHECK_THROWS(Product<>::value(dm02, sm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixSparseMatrixProductQuickTest<float> dense_matrix_sparse_matrix_product_quick_test_float("float");
DenseMatrixSparseMatrixProductQuickTest<double> dense_matrix_sparse_matrix_product_quick_test_double("double");

template <typename Tag_, typename DataType_>
class SparseMatrixDenseMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_matrix_product_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size + 3, size + 2, size / 8 + 1);
                DenseMatrix<DataType_> dm2(size + 2, size + 1);
                DenseMatrix<DataType_> dm3(size + 3 , size + 1, DataType_(0));
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 14 == 0)
                        *i = DataType_(i.index()) / 1.23456;
                    if (i.index() % 21 == 0)
                        *i = DataType_(i.index()) / 1.23456;
                }
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index()) / 3.9876;
                }

                typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
                for (unsigned int s(0) ; s < sm1.rows() ; ++s)
                {
                    const SparseVector<DataType_> sm1_row(sm1[s]);
                    for (unsigned int t(0); t < dm2.columns() ; ++t)
                    {
                        const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                        *i = DotProduct<>::value(sm1_row, dm2_column);
                        ++i;
                    }
                }
                DenseMatrix<DataType_> prod(Product<Tag_>::value(sm1, dm2));

                for (typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                        j(dm3.begin_elements()) ; i != i_end ; ++i,  ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 10 * *i * std::numeric_limits<DataType_>::epsilon());
                }
            }

            SparseMatrix<DataType_> sm01(3, 4, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(Product<Tag_>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductTest<tags::CPU, float> sparse_matrix_dense_matrix_product_test_float("float");
SparseMatrixDenseMatrixProductTest<tags::CPU, double> sparse_matrix_dense_matrix_product_test_double("double");
SparseMatrixDenseMatrixProductTest<tags::CPU::MultiCore, float> mc_sparse_matrix_dense_matrix_product_test_float("MC float");
SparseMatrixDenseMatrixProductTest<tags::CPU::MultiCore, double> mc_sparse_matrix_dense_matrix_product_test_double("MC double");
#ifdef HONEI_SSE
SparseMatrixDenseMatrixProductTest<tags::CPU::SSE, float> sse_sparse_matrix_dense_matrix_product_test_float("SSE float");
SparseMatrixDenseMatrixProductTest<tags::CPU::SSE, double> sse_sparse_matrix_dense_matrix_product_test_double("SSE double");
SparseMatrixDenseMatrixProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_sparse_matrix_dense_matrix_product_test_float("MC SSE float");
SparseMatrixDenseMatrixProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_sparse_matrix_dense_matrix_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
SparseMatrixDenseMatrixProductTest<tags::Cell, float> cell_sparse_matrix_dense_matrix_product_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_matrix_product_quick_test<" + type + ">")
    {
        register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size(11);
            SparseMatrix<DataType_> sm1(size + 3, size + 2, size / 8 + 1);
            DenseMatrix<DataType_> dm2(size + 2, size + 1);
            DenseMatrix<DataType_> dm3(size + 3 , size + 1, DataType_(0));
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 14 == 0)
                    *i = DataType_(i.index()) / 1.23456;
                if (i.index() % 21 == 0)
                    *i = DataType_(i.index()) / 1.23456;
            }
            for (typename DenseMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(i.index()) / 3.9876;
            }

            typename DenseMatrix<DataType_>::ElementIterator i(dm3.begin_elements());
            for (unsigned int s(0) ; s < sm1.rows() ; ++s)
            {
                const SparseVector<DataType_> sm1_row(sm1[s]);
                for (unsigned int t(0); t < dm2.columns() ; ++t)
                {
                    const DenseVectorSlice<DataType_> dm2_column(dm2.column(t));
                    *i = DotProduct<>::value(sm1_row, dm2_column);
                    ++i;
                }
            }
            DenseMatrix<DataType_> prod(Product<Tag_>::value(sm1, dm2));

            for (typename DenseMatrix<DataType_>::ConstElementIterator i(prod.begin_elements()), i_end(prod.end_elements()),
                    j(dm3.begin_elements()) ; i != i_end ; ++i,  ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 10 * *i * std::numeric_limits<DataType_>::epsilon());
            }

            SparseMatrix<DataType_> sm01(3, 4, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(Product<Tag_>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductQuickTest<tags::CPU, float> sparse_matrix_dense_matrix_product_quick_test_float("float");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU, double> sparse_matrix_dense_matrix_product_quick_test_double("double");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_dense_matrix_product_quick_test_float("MC float");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_dense_matrix_product_quick_test_double("MC double");
#ifdef HONEI_SSE
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::SSE, float> sse_sparse_matrix_dense_matrix_product_quick_test_float("SSE float");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::SSE, double> sse_sparse_matrix_dense_matrix_product_quick_test_double("SSE double");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_sparse_matrix_dense_matrix_product_quick_test_float("MC SSE float");
SparseMatrixDenseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_sparse_matrix_dense_matrix_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
SparseMatrixDenseMatrixProductQuickTest<tags::Cell, float> cell_sparse_matrix_dense_matrix_product_quick_test_float("Cell float");
#endif

template <typename DataType_>
class SparseMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_product_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                    sm2(size, size+1, size / 7 + 1), sm3(size + 1, size + 1, size / 7 + 1);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                        i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                        i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm3.begin_elements()),
                        i_end(sm3.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 6 * size;
                }

                SparseMatrix<DataType_> prod(Product<>::value(sm1, sm2));

                TEST_CHECK_EQUAL(prod, sm3);
            }

            SparseMatrix<DataType_> sm01(3, 3, 1), sm02(3, 4, 1);

            TEST_CHECK_THROWS(Product<>::value(sm02, sm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixProductTest<float> sparse_matrix_product_test_float("float");
SparseMatrixProductTest<double> sparse_matrix_product_test_double("double");

template <typename DataType_>
class SparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_product_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(2);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                sm2(size, size+1, size / 7 + 1), sm3(size + 1, size + 1, size / 7 + 1);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                    i_end(sm2.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm3.begin_elements()),
                    i_end(sm3.end_elements()) ; i != i_end ; ++i)
            {
                *i = 6 * size;
            }

            SparseMatrix<DataType_> prod(Product<>::value(sm1, sm2));

            TEST_CHECK_EQUAL(prod, sm3);

            SparseMatrix<DataType_> sm01(3, 3, 1), sm02(3, 4, 1);

            TEST_CHECK_THROWS(Product<>::value(sm02, sm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixProductQuickTest<float> sparse_matrix_product_quick_test_float("float");
SparseMatrixProductQuickTest<double> sparse_matrix_product_quick_test_double("double");

template <typename DataType_>
class BandedMatrixDenseMatrixProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_product_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                DenseVector<DataType_> dv4(size, DataType_(6));
                DenseMatrix<DataType_> dm1(size, size, DataType_(0));
                for(unsigned long i=0 ; i < size ; ++i)
                {
                    dm1[i][i] = DataType_(3);
                }
                BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);
                bm1.insert_band(1, dv3.copy());
                bm1.insert_band(-2, dv4.copy());

                DenseVector<DataType_> dv5(size, DataType_(18));
                DenseVector<DataType_> dv6(size, DataType_(18));

                bm3.insert_band(1, dv5);
                bm3.insert_band(-2, dv6);
                DenseMatrix<DataType_> prod(Product<>::value(bm1, dm1));

                typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
                for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
                }

                BandedMatrix<DataType_> bm01(5);
                DenseMatrix<DataType_> dm02(6, 5);

                TEST_CHECK_THROWS(Product<>::value(bm01, dm02), MatrixRowsDoNotMatch);
            }
        }
};
BandedMatrixDenseMatrixProductTest<float> banded_matrix_dense_matrix_product_test_float("float");
BandedMatrixDenseMatrixProductTest<double> banded_matrix_dense_matrix_product_test_double("double");

template <typename DataType_>
class BandedMatrixDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_product_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            DenseVector<DataType_> dv4(size, DataType_(6));
            DenseMatrix<DataType_> dm1(size, size, DataType_(0));
            for(unsigned long i=0 ; i < size ; ++i)
            {
                dm1[i][i] = DataType_(3);
            }
            BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);
            bm1.insert_band(1, dv3.copy());
            bm1.insert_band(-2, dv4.copy());

            DenseVector<DataType_> dv5(size, DataType_(18));
            DenseVector<DataType_> dv6(size, DataType_(18));

            bm3.insert_band(1, dv5);
            bm3.insert_band(-2, dv6);
            DenseMatrix<DataType_> prod(Product<>::value(bm1, dm1));

            typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
            for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
            }

            BandedMatrix<DataType_> bm01(5);
            DenseMatrix<DataType_> dm02(6, 5);

            TEST_CHECK_THROWS(Product<>::value(bm01, dm02), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseMatrixProductQuickTest<float> banded_matrix_dense_matrix_product_quick_test_float("float");
BandedMatrixDenseMatrixProductQuickTest<double> banded_matrix_dense_matrix_product_quick_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixProductTest :
    public BaseTest
{
    public:
        BandedMatrixSparseMatrixProductTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                DenseVector<DataType_> dv4(size, DataType_(6));

                SparseMatrix<DataType_> sm1(size, size);
                sm1[0][size-1] = DataType_(3);
                sm1[size-1][0] = DataType_(3);

                for(unsigned long i=0; i < size ; ++i)
                {
                    sm1[i][i] = DataType_(3);
                }

                BandedMatrix<DataType_> bm1(size, dv1);
                bm1.insert_band(1, dv3.copy());
                bm1.insert_band(-2, dv4.copy());

                // Create a DenseMatrix that equals the BandedMatrix:

                DenseMatrix<DataType_> dm1(size, size, DataType_(0));
                for (unsigned long i = 0; i < size; ++i)
                {
                    for (unsigned long j = 0; j < size; ++j)
                    {
                        if(i == j)
                            dm1[i][i] = DataType_(2);

                        if (j == i+1)
                            dm1[i][j] = DataType_(6);

                        if (j == i-2)
                            dm1[i][j] = DataType_(6);
                    }
                }

                //Calculate products for dense * sparse and banded * sparse and check equality.
                DenseMatrix<DataType_> prod(Product<>::value(bm1, sm1));
                DenseMatrix<DataType_> prod2(Product<>::value(dm1, sm1));
                TEST_CHECK_EQUAL(prod, prod2);

                BandedMatrix<DataType_> bm01(5);
                SparseMatrix<DataType_> sm02(6, 5);

                TEST_CHECK_THROWS(Product<>::value(bm01, sm02), MatrixRowsDoNotMatch);
            }
        }
};
BandedMatrixSparseMatrixProductTest<float> banded_matrix_sparse_matrix_product_test_float("float");
BandedMatrixSparseMatrixProductTest<double> banded_matrix_sparse_matrix_product_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            DenseVector<DataType_> dv4(size, DataType_(6));

            SparseMatrix<DataType_> sm1(size, size);
            sm1[0][size-1] = DataType_(3);
            sm1[size-1][0] = DataType_(3);

            for(unsigned long i=0; i < size ; ++i)
            {
                sm1[i][i] = DataType_(3);
            }

            BandedMatrix<DataType_> bm1(size, dv1);
            bm1.insert_band(1, dv3.copy());
            bm1.insert_band(-2, dv4.copy());

            // Create a DenseMatrix that equals the BandedMatrix:

            DenseMatrix<DataType_> dm1(size, size, DataType_(0));
            for (unsigned long i = 0; i < size; ++i)
            {
                for (unsigned long j = 0; j < size; ++j)
                {
                    if(i == j)
                        dm1[i][i] = DataType_(2);

                    if (j == i+1)
                        dm1[i][j] = DataType_(6);

                    if (j == i-2)
                        dm1[i][j] = DataType_(6);
                }
            }

            //Calculate products for dense * sparse and banded * sparse and check equality.
            DenseMatrix<DataType_> prod(Product<>::value(bm1, sm1));
            DenseMatrix<DataType_> prod2(Product<>::value(dm1, sm1));
            TEST_CHECK_EQUAL(prod, prod2);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5);

            TEST_CHECK_THROWS(Product<>::value(bm01, sm02), MatrixRowsDoNotMatch);
        }
};
BandedMatrixSparseMatrixProductQuickTest<float> banded_matrix_sparse_matrix_product_quick_test_float("float");
BandedMatrixSparseMatrixProductQuickTest<double> banded_matrix_sparse_matrix_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixBandedMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixBandedMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_banded_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size, DataType_(0));
                for (unsigned long i = 0; i < size; ++i)
                {
                    for (unsigned long j = 0; j < size; ++j)
                    {
                        if(i == j)
                            dm1[i][i] = DataType_(2);

                        if (j == i+1)
                            dm1[i][j] = DataType_(6);

                        if (j == i-2)
                            dm1[i][j] = DataType_(6);
                    }
                }
                DenseVector<DataType_> dv1 (size, DataType_(3));
                DenseVector<DataType_> dv3 (size, DataType_(6));

                BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);

                DenseVector<DataType_> dv5 (size, DataType_(18));

                bm3.insert_band(1, dv5);
                bm3.insert_band(-2, dv5);
                DenseMatrix<DataType_> prod(Product<>::value(dm1, bm1));

                typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
                for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
                }


                DenseMatrix<DataType_> dm03(5, 6);
                BandedMatrix<DataType_> bm01(5);

                TEST_CHECK_THROWS(Product<>::value(dm03, bm01), MatrixRowsDoNotMatch);
            }
        }
};
DenseMatrixBandedMatrixProductTest<float> dense_matrix_banded_matrix_product_test_float("float");
DenseMatrixBandedMatrixProductTest<double> dense_matrix_banded_matrix_product_test_double("double");

template <typename DataType_>
class DenseMatrixBandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixBandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_banded_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);

            DenseMatrix<DataType_> dm1(size, size, DataType_(0));
            for (unsigned long i = 0; i < size; ++i)
            {
                for (unsigned long j = 0; j < size; ++j)
                {
                    if(i == j)
                        dm1[i][i] = DataType_(2);

                    if (j == i+1)
                        dm1[i][j] = DataType_(6);

                    if (j == i-2)
                        dm1[i][j] = DataType_(6);
                }
            }
            DenseVector<DataType_> dv1 (size, DataType_(3));
            DenseVector<DataType_> dv3 (size, DataType_(6));

            BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);

            DenseVector<DataType_> dv5 (size, DataType_(18));

            bm3.insert_band(1, dv5);
            bm3.insert_band(-2, dv5);
            DenseMatrix<DataType_> prod(Product<>::value(dm1, bm1));

            typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
            for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
            }


            DenseMatrix<DataType_> dm03(5, 6);
            BandedMatrix<DataType_> bm01(5);

            TEST_CHECK_THROWS(Product<>::value(dm03, bm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixBandedMatrixProductQuickTest<float> dense_matrix_banded_matrix_product_quick_test_float("float");
DenseMatrixBandedMatrixProductQuickTest<double> dense_matrix_banded_matrix_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixBandedMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixBandedMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_banded_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size);
                for (unsigned long i = 0; i < size; ++i)
                {
                    for (unsigned long j = 0; j < size; ++j)
                    {
                        if(i == j)
                            sm1[i][i] = DataType_(2);

                        if (j == i+1)
                            sm1[i][j] = DataType_(6);

                        if (j == i-2)
                            sm1[i][j] = DataType_(6);
                    }
                }
                DenseVector<DataType_> dv1(size, DataType_(3));
                DenseVector<DataType_> dv3(size, DataType_(6));

                BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);

                DenseVector<DataType_> dv5(size, DataType_(18));
                DenseVector<DataType_> dv6(size, DataType_(18));
                bm3.insert_band(1, dv5.copy());
                bm3.insert_band(-2, dv5.copy());
                DenseMatrix<DataType_> prod(Product<>::value(sm1, bm1));

                typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
                for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
                }

                SparseMatrix<DataType_> sm03(5, 6);
                BandedMatrix<DataType_> bm01(5);

                TEST_CHECK_THROWS(Product<>::value(sm03, bm01), MatrixRowsDoNotMatch);
            }
        }
};
SparseMatrixBandedMatrixProductTest<float> sparse_matrix_banded_matrix_product_test_float("float");
SparseMatrixBandedMatrixProductTest<double> sparse_matrix_banded_matrix_product_test_double("double");

template <typename DataType_>
class SparseMatrixBandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixBandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_banded_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);

            SparseMatrix<DataType_> sm1(size, size);
            for (unsigned long i = 0; i < size; ++i)
            {
                for (unsigned long j = 0; j < size; ++j)
                {
                    if(i == j)
                        sm1[i][i] = DataType_(2);

                    if (j == i+1)
                        sm1[i][j] = DataType_(6);

                    if (j == i-2)
                        sm1[i][j] = DataType_(6);
                }
            }
            DenseVector<DataType_> dv1(size, DataType_(3));
            DenseVector<DataType_> dv3(size, DataType_(6));

            BandedMatrix<DataType_> bm1(size, dv1), bm3(size, dv3);

            DenseVector<DataType_> dv5(size, DataType_(18));
            DenseVector<DataType_> dv6(size, DataType_(18));
            bm3.insert_band(1, dv5.copy());
            bm3.insert_band(-2, dv5.copy());
            DenseMatrix<DataType_> prod(Product<>::value(sm1, bm1));

            typename BandedMatrix<DataType_>::ConstElementIterator k(bm3.begin_elements());
            for (typename DenseMatrix<DataType_>::ConstElementIterator l(prod.begin_elements()), l_end(prod.end_elements()) ; l != l_end ; ++l, ++k)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*k, *l, std::numeric_limits<DataType_>::epsilon());
            }

            SparseMatrix<DataType_> sm03(5, 6);
            BandedMatrix<DataType_> bm01(5);

            TEST_CHECK_THROWS(Product<>::value(sm03, bm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixBandedMatrixProductQuickTest<float> sparse_matrix_banded_matrix_product_quick_test_float("float");
SparseMatrixBandedMatrixProductQuickTest<double> sparse_matrix_banded_matrix_product_quick_test_double("double");

