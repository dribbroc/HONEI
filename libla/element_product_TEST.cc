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

#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/element_product.hh>
#include <libla/vector_error.cc>
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
                DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                    dv3(size, DataType_(6));
                DenseVector<DataType_> prod(ElementProduct<>::value(dv1, dv2));

                TEST_CHECK_EQUAL(prod, dv3);
            }

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<>::value(dv02, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorElementProductTest<tags::CPU, float> dense_vector_elementwise_product_test_float("float");
DenseVectorElementProductTest<tags::CPU, double> dense_vector_elementwise_product_test_double("double");
#ifdef HONEI_SSE
DenseVectorElementProductTest<tags::CPU::SSE, float> sse_dense_vector_elementwise_product_test_float("sse float");
DenseVectorElementProductTest<tags::CPU::SSE, double> sse_dense_vector_elementwise_product_test_double("sse double");
#endif
#ifdef HONEI_CELL
//DenseVectorElementProductTest<tags::Cell, float> cell_dense_vector_element_product_test_float("Cell float");
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
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                dv3(size, DataType_(6));
            DenseVector<DataType_> prod(ElementProduct<Tag_>::value(dv1, dv2));

            TEST_CHECK_EQUAL(prod, dv3);


            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<Tag_>::value(dv02, dv01), VectorSizeDoesNotMatch);

        }
};
DenseVectorElementProductQuickTest<tags::CPU, float> dense_vector_elementwise_product_quick_test_float("float");
DenseVectorElementProductQuickTest<tags::CPU, double> dense_vector_elementwise_product_quick_test_double("double");
#ifdef HONEI_SSE
DenseVectorElementProductQuickTest<tags::CPU, float> sse_dense_vector_elementwise_product_quick_test_float("sse float");
DenseVectorElementProductQuickTest<tags::CPU, double> sse_dense_vector_elementwise_product_quick_test_double("sse double");
#endif
#ifdef HONEI_CELL
DenseVectorElementProductQuickTest<tags::Cell, float> cell_dense_vector_element_product_quick_test_float("Cell float");
#endif


template <typename DataType_>
class SparseVectorDenseVectorElementProductTest :
    public BaseTest
{
    public:
        SparseVectorDenseVectorElementProductTest(const std::string & type) :
            BaseTest("sparse_vector_dense_vector_elementwise_product_test<" + type + ">")
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

                SparseVector<DataType_> prod(ElementProduct<>::value(sv1, dv2));

                TEST_CHECK_EQUAL(prod, sv3);
            }

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(ElementProduct<>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementProductTest<float> sparse_vector_dense_vector_elementwise_product_test_float("float");
SparseVectorDenseVectorElementProductTest<double> sparse_vector_dense_vector_elementwise_product_test_double("double");

template <typename DataType_>
class SparseVectorDenseVectorElementProductQuickTest :
    public QuickTest
{
    public:
        SparseVectorDenseVectorElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_vector_dense_vector_elementwise_product_quick_test<" + type + ">")
        {
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
            SparseVector<DataType_> prod(ElementProduct<>::value(sv1, dv2));

            TEST_CHECK_EQUAL(prod, sv3);

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(ElementProduct<>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementProductQuickTest<float> sparse_vector_dense_vector_elementwise_product_quick_test_float("float");
SparseVectorDenseVectorElementProductQuickTest<double> sparse_vector_dense_vector_elementwise_product_quick_test_double("double");

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
                SparseVector<DataType_> prod(ElementProduct<>::value(sv1, sv2));

                TEST_CHECK_EQUAL(prod, sv3);
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
            SparseVector<DataType_> prod(ElementProduct<>::value(sv1, sv2));

            TEST_CHECK_EQUAL(prod, sv3);

            SparseVector<DataType_> sv01(3, 1), sv02(4, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sv02, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorElementProductQuickTest<float> sparse_vector_elementwise_product_quick_test_float("float");
SparseVectorElementProductQuickTest<double> sparse_vector_elementwise_product_quick_test_double("double");
template <typename DataType_>
class BandedMatrixDenseMatrixElementProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixElementProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size, DataType_(3));
                DenseVector<DataType_> dv2(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3);

                BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm2, dm1));

                TEST_CHECK_EQUAL(prod, bm3);
            }

            DenseMatrix<DataType_> dm01(5, 6), dm03(6, 5);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, dm01), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixDenseMatrixElementProductTest<float> banded_matrix_dense_matrix_elementwise_product_test_float("float");
BandedMatrixDenseMatrixElementProductTest<double> banded_matrix_dense_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixDenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size, DataType_(3));
            DenseVector<DataType_> dv2(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3);

            BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm2, dm1));

            TEST_CHECK_EQUAL(prod, bm3);

            DenseMatrix<DataType_> dm01(6, 5), dm03(5, 6);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, dm01), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixDenseMatrixElementProductQuickTest<float> banded_matrix_dense_matrix_elementwise_product_quick_test_float("float");
BandedMatrixDenseMatrixElementProductQuickTest<double> banded_matrix_dense_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class BandedMatrixElementProductTest :
    public BaseTest
{
    public:
        BandedMatrixElementProductTest(const std::string & type) :
            BaseTest("banded_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));
                DenseVector<DataType_> dv2(size, DataType_(2));
                DenseVector<DataType_> dv3(size, DataType_(6));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);
                BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm1, bm2));

                TEST_CHECK_EQUAL(prod, bm3);
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixElementProductTest<float> banded_matrix_elementwise_product_test_float("float");
BandedMatrixElementProductTest<double> banded_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> dv1(size, DataType_(3));
            DenseVector<DataType_> dv2(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(6));
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);
            BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm1, bm2));

            TEST_CHECK_EQUAL(prod, bm3);

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixElementProductQuickTest<float> banded_matrix_elementwise_product_quick_test_float("float");
BandedMatrixElementProductQuickTest<double> banded_matrix_elementwise_product_quick_test_double("double");

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

                BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm1, sm2));

                TEST_CHECK_EQUAL(prod, bm3);
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

            BandedMatrix<DataType_> & prod(ElementProduct<>::value(bm1, sm2));

            TEST_CHECK_EQUAL(prod, bm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm02), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(bm01, sm03), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixSparseMatrixElementProductQuickTest<float> banded_matrix_sparse_matrix_elementwise_product_quick_test_float("float");
BandedMatrixSparseMatrixElementProductQuickTest<double> banded_matrix_sparse_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixElementProductTest :
    public BaseTest
{
    public:
        DenseMatrixElementProductTest(const std::string & type) :
            BaseTest("dense_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(3)),
                    dm3(size+1, size, DataType_(6));
                DenseMatrix<DataType_> & prod(ElementProduct<>::value(dm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            DenseMatrix<DataType_> dm01(3, 2, static_cast<DataType_>(1)), dm02(4, 3, static_cast<DataType_>(1)),
                dm03(4, 2, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(ElementProduct<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementProductTest<float> dense_matrix_elementwise_product_test_float("float");
DenseMatrixElementProductTest<double> dense_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class DenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(3)),
                dm3(size+1, size, DataType_(6));
            DenseMatrix<DataType_> & prod(ElementProduct<>::value(dm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            DenseMatrix<DataType_> dm01(3, 2, static_cast<DataType_>(1)), dm02(4, 3, DataType_(1)),
                dm03(4, 2, DataType_(1));

            TEST_CHECK_THROWS(ElementProduct<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementProductQuickTest<float> dense_matrix_elementwise_product_quick_test_float("float");
DenseMatrixElementProductQuickTest<double> dense_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixDenseMatrixElementProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseMatrixElementProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
                SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), j(sm2.begin_elements()),
                    k(sm3.begin_elements()) ; j != j_end ; ++j, ++k)
                {
                    if (j.index() % 7 == 0)
                    {
                        *j = DataType_(3);
                        *k = DataType_(6);
                    }
                }

                SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm2, dm1));

                TEST_CHECK_EQUAL(prod, sm3);
            }

            DenseMatrix<DataType_> dm03(6, 5);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm01, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm02, dm03), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDenseMatrixElementProductTest<float> sparse_matrix_dense_matrix_elementwise_product_test_float("float");
SparseMatrixDenseMatrixElementProductTest<double> sparse_matrix_dense_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class SparseMatrixDenseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2));
            SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), j(sm2.begin_elements()),
                k(sm3.begin_elements()) ; j != j_end ; ++j, ++k)
            {
                if (j.index() % 7 == 0)
                {
                    *j = DataType_(3);
                    *k = DataType_(6);
                }
            }

            SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm2, dm1));

            TEST_CHECK_EQUAL(prod, sm3);

            DenseMatrix<DataType_> dm03(6, 5);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm01, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm02, dm03), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDenseMatrixElementProductQuickTest<float> sparse_matrix_dense_matrix_elementwise_product_quick_test_float("float");
SparseMatrixDenseMatrixElementProductQuickTest<double> sparse_matrix_dense_matrix_elementwise_product_quick_test_double("double");

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

                SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm1, bm2));

                TEST_CHECK_EQUAL(prod, sm3);
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

            SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm1, bm2));

            TEST_CHECK_EQUAL(prod, sm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, bm01), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixBandedMatrixElementProductQuickTest<float> sparse_matrix_banded_matrix_elementwise_product_quick_test_float("float");
SparseMatrixBandedMatrixElementProductQuickTest<double> sparse_matrix_banded_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixElementProductTest :
    public BaseTest
{
    public:
        SparseMatrixElementProductTest(const std::string & type) :
            BaseTest("sparse_matrix_elementwise_product_test<" + type + ">")
        {
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

                SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm1, sm2));

                TEST_CHECK_EQUAL(prod, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixElementProductTest<float> sparse_matrix_elementwise_product_test_float("float");
SparseMatrixElementProductTest<double> sparse_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class SparseMatrixElementProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
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

            SparseMatrix<DataType_> & prod(ElementProduct<>::value(sm1, sm2));

            TEST_CHECK_EQUAL(prod, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(ElementProduct<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixElementProductQuickTest<float> sparse_matrix_elementwise_product_quick_test_float("float");
SparseMatrixElementProductQuickTest<double> sparse_matrix_elementwise_product_quick_test_double("double");
