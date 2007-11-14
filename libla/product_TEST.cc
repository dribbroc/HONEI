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

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/product.hh>
#include <libla/matrix_error.cc>
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
                DenseVector<DataType_> dv1(size, DataType_(2));
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseVector<DataType_> dv4(size, DataType_(2));
                DenseVector<DataType_> dv5(dv4.copy());
                bm1.insert_band(1, dv4);
                bm1.insert_band(-1, dv5);

                DenseVector<DataType_> dv2(size, DataType_(3)), dv3(size, DataType_(6 * 3));
                (dv3)[0]= DataType_(12);
                (dv3)[size-1]= DataType_(12);
                DenseVector<DataType_> prod (Product<Tag_>::value(bm1, dv2));

                TEST_CHECK_EQUAL(prod, dv3);
            }

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixDenseVectorProductTest<tags::CPU, float> banded_matrix_dense_vector_product_test_float("float");
BandedMatrixDenseVectorProductTest<tags::CPU, double> banded_matrix_dense_vector_product_test_double("double");
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductTest<tags::CPU::SSE, float> sse_banded_matrix_dense_vector_product_test_float("SSE float");
BandedMatrixDenseVectorProductTest<tags::CPU::SSE, double> sse_banded_matrix_dense_vector_product_test_double("SSE double");
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
            unsigned long size(43);
            DenseVector<DataType_> dv1(size, DataType_(2));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseVector<DataType_> dv4(size, DataType_(2));
            DenseVector<DataType_> dv5(dv4.copy());
            bm1.insert_band(1, dv4);
            bm1.insert_band(-1, dv5);

            DenseVector<DataType_> dv2(size, DataType_(3)), dv3(size, DataType_(6 * 3));
            (dv3)[0]= DataType_(12);
            (dv3)[size-1]= DataType_(12);
            DenseVector<DataType_> prod(Product<Tag_>::value(bm1, dv2));

            TEST_CHECK_EQUAL(prod, dv3);

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
BandedMatrixDenseVectorProductQuickTest<tags::CPU, float> banded_matrix_dense_vector_product_quick_test_float("float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU, double> banded_matrix_dense_vector_product_quick_test_double("double");
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductQuickTest<tags::CPU::SSE, float> sse_banded_matrix_dense_vector_product_quick_test_float("SSE float");
BandedMatrixDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_banded_matrix_dense_vector_product_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductQuickTest<tags::Cell, float> cell_banded_matrix_dense_vector_product_quick_test_float("CELL float");
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
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
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
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 5 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
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
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
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
DenseMatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_matrix_dense_vector_product_test_float("MC SSE float");
DenseMatrixDenseVectorProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_matrix_dense_vector_product_test_double("MC SSE double");
DenseMatrixDenseVectorProductTest<tags::CPU::SSE, float> sse_dense_matrix_dense_vector_product_test_float("SSE float");
DenseMatrixDenseVectorProductTest<tags::CPU::SSE, double> sse_dense_matrix_dense_vector_product_test_double("SSE double");
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
DenseMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_matrix_dense_vector_product_quick_test_float("MC SSE float");
DenseMatrixDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_matrix_dense_vector_product_quick_test_double("MC SSE double");
DenseMatrixDenseVectorProductQuickTest<tags::CPU::SSE, float> sse_dense_matrix_dense_vector_product_quick_test_float("SSE float");
DenseMatrixDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_dense_matrix_dense_vector_product_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixDenseVectorProductQuickTest<tags::Cell, float> cell_dense_matrix_dense_vector_product_quick_test_float("Cell float");
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
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(6 * (size / 10 + 1));
                }
                DenseVector<DataType_> prod(Product<>::value(dm1, sv1));

                TEST_CHECK_EQUAL(prod, sv2);
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
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(12);
            }
            DenseVector<DataType_> prod(Product<>::value(dm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

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
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                SparseVector<DataType_> prod(Product<>::value(sm1, dv1));

                TEST_CHECK_EQUAL(prod, dv2);
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
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
            SparseVector<DataType_> prod(Product<>::value(sm1, dv1));

            TEST_CHECK_EQUAL(prod, dv2);

            SparseMatrix<DataType_> sm01(4, 3, 1);
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(Product<>::value(sm01, dv01), VectorSizeDoesNotMatch);

        }
};
SparseMatrixDenseVectorProductQuickTest<float> sparse_matrix_dense_vector_product_quick_test_float("float");
SparseMatrixDenseVectorProductQuickTest<double> sparse_matrix_dense_vector_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        SparseMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("sparse__matrix_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
                {
                    *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
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
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
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

            for(int i=0; i < size ; ++i)
            {
                sm1[i][i] = DataType_(3);
            }

            //Calculate products for banded * banded and banded * sparse and check equality.
            BandedMatrix<DataType_> prod(Product<>::value(bm1, bm2));
            DenseMatrix<DataType_> prod2(Product<>::value(bm1, sm1));
            TEST_CHECK_EQUAL(prod, prod2);

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Product<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixProductQuickTest<float> banded_matrix_product_quick_test_float("float");
BandedMatrixProductQuickTest<double> banded_matrix_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size, size+1, DataType_(3)),
                    dm3(size+1, size+1, DataType_(6 * size));
                DenseMatrix<DataType_> prod(Product<>::value(dm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1)), dm02(3, 3, DataType_(1));

            TEST_CHECK_THROWS(Product<>::value(dm01, dm02), MatrixRowsDoNotMatch);
        }
};
DenseMatrixProductTest<float> dense_matrix_product_test_float("float");
DenseMatrixProductTest<double> dense_matrix_product_test_double("double");

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
            unsigned long size(31);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size, size+1, DataType_(3)),
                dm3(size+1, size+1, DataType_(6 * size));

            DenseMatrix<DataType_> prod(Product<Tag_>::value(dm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(3, 4, DataType_(1));

            TEST_CHECK_THROWS(Product<Tag_>::value(dm02, dm01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixProductQuickTest<tags::CPU, float> dense_matrix_product_quick_test_float("float");
DenseMatrixProductQuickTest<tags::CPU, double> dense_matrix_product_quick_test_double("double");
#ifdef HONEI_CELL
DenseMatrixProductQuickTest<tags::Cell, float> cell_dense_matrix_product_quick_test_float("Cell float");
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
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
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
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
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

template <typename DataType_>
class SparseMatrixDenseMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm2(size, size+1, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }
                DenseMatrix<DataType_> prod(Product<>::value(sm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            SparseMatrix<DataType_> sm01(3, 4, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(Product<>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductTest<float> sparse_matrix_dense_matrix_product_test_float("float");
SparseMatrixDenseMatrixProductTest<double> sparse_matrix_dense_matrix_product_test_double("double");

template <typename DataType_>
class SparseMatrixDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseMatrix<DataType_> dm2(size, size+1, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }
            DenseMatrix<DataType_> prod(Product<>::value(sm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            SparseMatrix<DataType_> sm01(3, 4, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(Product<>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductQuickTest<float> sparse_matrix_dense_matrix_product_quick_test_float("float");
SparseMatrixDenseMatrixProductQuickTest<double> sparse_matrix_dense_matrix_product_quick_test_double("double");

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
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                    i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm3.begin_elements()),
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
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                i_end(sm2.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm3.begin_elements()),
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
                for(int i=0 ; i < size ; ++i)
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

                TEST_CHECK_EQUAL(prod, bm3);

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
            for(int i=0 ; i < size ; ++i)
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

            TEST_CHECK_EQUAL(prod, bm3);

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

                for(int i=0; i < size ; ++i)
                {
                    sm1[i][i] = DataType_(3);
                }

                BandedMatrix<DataType_> bm1(size, dv1);
                bm1.insert_band(1, dv3.copy());
                bm1.insert_band(-2, dv4.copy());

                // Create a DenseMatrix that equals the BandedMatrix:

                DenseMatrix<DataType_> dm1(size, size, DataType_(0));
                for (int i = 0; i < size; ++i)
                {
                    for (int j = 0; j < size; ++j)
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

            for(int i=0; i < size ; ++i)
            {
                sm1[i][i] = DataType_(3);
            }

            BandedMatrix<DataType_> bm1(size, dv1);
            bm1.insert_band(1, dv3.copy());
            bm1.insert_band(-2, dv4.copy());

            // Create a DenseMatrix that equals the BandedMatrix:

            DenseMatrix<DataType_> dm1(size, size, DataType_(0));
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < size; ++j)
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
                for (int i = 0; i < size; ++i)
                {
                    for (int j = 0; j < size; ++j)
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

                TEST_CHECK_EQUAL(prod, bm3);

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
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < size; ++j)
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

            TEST_CHECK_EQUAL(prod, bm3);

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
                for (int i = 0; i < size; ++i)
                {
                    for (int j = 0; j < size; ++j)
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

                TEST_CHECK_EQUAL(prod, bm3);

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
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < size; ++j)
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

            TEST_CHECK_EQUAL(prod, bm3);

            SparseMatrix<DataType_> sm03(5, 6);
            BandedMatrix<DataType_> bm01(5);

            TEST_CHECK_THROWS(Product<>::value(sm03, bm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixBandedMatrixProductQuickTest<float> sparse_matrix_banded_matrix_product_quick_test_float("float");
SparseMatrixBandedMatrixProductQuickTest<double> sparse_matrix_banded_matrix_product_quick_test_double("double");
