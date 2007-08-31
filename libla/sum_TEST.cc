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
#include <libla/matrix_error.hh>
#include <libla/norm.hh>
#include <libla/sum.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;

// Test cases for matrix operations

template <typename DataType_>
class BandedMatrixDenseMatrixSumTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1(new DenseVector<DataType_>(size, DataType_(2)));
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseMatrix<DataType_> dm2(size, size, DataType_(1)), dm3(size, size, DataType_(1));

                typename MutableMatrix<DataType_>::ElementIterator k(dm3.begin_elements());
                for (typename Matrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                    i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
                {
                    if (*i != DataType_(0))
                    {
                        *k = DataType_(3);
                    }
                }
                DenseMatrix<DataType_> & sum(Sum<>::value(dm2, bm1));

                TEST_CHECK_EQUAL(sum, dm3);
            }

            BandedMatrix<DataType_> bm01(5);
            DenseMatrix<DataType_> dm02(6, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixDenseMatrixSumTest<float> banded_matrix_dense_matrix_sum_test_float("float");
BandedMatrixDenseMatrixSumTest<double> banded_matrix_dense_matrix_sum_test_double("double");

template <typename DataType_>
class BandedMatrixDenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> * dv1(new DenseVector<DataType_>(size, DataType_(2)));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseMatrix<DataType_> dm2(size, size, DataType_(1)), dm3(size, size, DataType_(1));

            typename MutableMatrix<DataType_>::ElementIterator k(dm3.begin_elements());
            for (typename Matrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
            {
                if (*i != DataType_(0))
                {
                    *k = DataType_(3);
                }
            }
            DenseMatrix<DataType_> & sum(Sum<>::value(dm2, bm1));

            TEST_CHECK_EQUAL(sum, dm3);

            BandedMatrix<DataType_> bm01(5);
            DenseMatrix<DataType_> dm02(6, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm02, bm01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseMatrixSumQuickTest<float> banded_matrix_dense_matrix_sum_quick_test_float("float");
BandedMatrixDenseMatrixSumQuickTest<double> banded_matrix_dense_matrix_sum_quick_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixSumTest :
    public BaseTest
{
    public:
        BandedMatrixSparseMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1(new DenseVector<DataType_>(size, DataType_(2)));
                BandedMatrix<DataType_> bm1(size, dv1);
                SparseMatrix<DataType_> sm2(size, size, size / 8 + 1),
                        sm3(size, size, size / 8 + 1);

                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                    i_end(sm2.end_elements()), k(sm3.begin_elements()) ;
                    i != i_end ; ++i, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(1);
                        *k = DataType_(1);
                    }
                }

                typename MutableMatrix<DataType_>::ElementIterator k(sm3.begin_elements());
                for (typename Matrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                    i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
                {
                    if (*i != DataType_(0))
                    {
                        *k = DataType_(3);
                    }
                }
                SparseMatrix<DataType_> & sum(Sum<>::value(sm2, bm1));

                TEST_CHECK_EQUAL(sum, sm3);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 5, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Sum<>::value(sm03, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(sm02, bm01), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixSparseMatrixSumTest<float> banded_matrix_sparse_matrix_sum_test_float("float");
BandedMatrixSparseMatrixSumTest<double> banded_matrix_sparse_matrix_sum_test_double("double");

template <typename DataType_>
class BandedMatrixSumTest :
    public BaseTest
{
    public:
        BandedMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3)));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
                BandedMatrix<DataType_> & sum(Sum<>::value(bm1, bm2));
                for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(sum.begin_bands()), ce_end(sum.end_bands()) ;
                    ce != ce_end ; ++ce)
                {
                    DenseVector<DataType_>  dv3 = *ce;
                    if (ce.index() != size - 1 )
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()),
                            i_end(dv3.end_elements()) ; i != i_end ; ++i)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                        }
                    }
                    else
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()),
                            i_end(dv3.end_elements()) ; i != i_end ; ++i)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                        }
                    }
                }
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Sum<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixSumTest<float> banded_matrix_sum_test_float("float");
BandedMatrixSumTest<double> banded_matrix_sum_test_double("double");

template <typename DataType_>
class BandedMatrixSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3)));
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
            BandedMatrix<DataType_> & sum(Sum<>::value(bm1, bm2));
            for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(sum.begin_bands()), ce_end(sum.end_bands()) ;
                ce != ce_end ; ++ce)
            {
                DenseVector<DataType_>  dv3 = *ce;
                if (ce.index() != size - 1 )
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()),
                        i_end(dv3.end_elements()) ; i != i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    }
                }
                else
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()),
                        i_end(dv3.end_elements()) ; i != i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                    }
                }
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Sum<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixSumQuickTest<float> banded_matrix_sum_quick_test_float("float");
BandedMatrixSumQuickTest<double> banded_matrix_sum_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixSumTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixSumTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(0)),
                        dm3(size, size + 1, DataType_(0));
                SparseMatrix<DataType_> sm2(size, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()),
                    i_end(dm1.end_elements()), j(sm2.begin_elements()), k(dm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_((i.index() +1)  / 1.23456789);
                        *k = DataType_((i.index() +1)  / 1.23456789);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_((i.index() +1) / 1.23456789);
                        *k = DataType_((i.index() +1) / 1.23456789);
                    }
                    if (i.index() % 10 == 0 && i.index() % 7 == 0)
                    {
                        *k = DataType_((i.index() +1) * 2 / 1.23456789);
                    }
                }
                DenseMatrix<DataType_> & sum(Sum<>::value(dm1, sm2));

                TEST_CHECK_EQUAL(sum, dm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DataType_> dm03(5, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixSumTest<float> dense_matrix_sparse_matrix_sum_test_float("float");
DenseMatrixSparseMatrixSumTest<double> dense_matrix_sparse_matrix_sum_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(0)),
                    dm3(size, size + 1, DataType_(0));
            SparseMatrix<DataType_> sm2(size, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()),
                i_end(dm1.end_elements()), j(sm2.begin_elements()), k(dm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_((i.index() +1)  / 1.23456789);
                    *k = DataType_((i.index() +1)  / 1.23456789);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_((i.index() +1) / 1.23456789);
                    *k = DataType_((i.index() +1) / 1.23456789);
                }
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                {
                    *k = DataType_((i.index() +1) * 2 / 1.23456789);
                }
            }
            DenseMatrix<DataType_> & sum(Sum<>::value(dm1, sm2));

            TEST_CHECK_EQUAL(sum, dm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DataType_> dm03(5, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixSumQuickTest<float> dense_matrix_sparse_matrix_sum_quick_test_float("float");
DenseMatrixSparseMatrixSumQuickTest<double> dense_matrix_sparse_matrix_sum_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSumTest :
    public BaseTest
{
    public:
        DenseMatrixSumTest(const std::string & type) :
            BaseTest("dense_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(-3)),
                    dm3(size, size + 1, DataType_(-1));
                DenseMatrix<DataType_> & sum(Sum<>::value(dm1, dm2));

                TEST_CHECK_EQUAL(sum, dm3);
            }

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSumTest<float> dense_matrix_sum_test_float("float");
DenseMatrixSumTest<double> dense_matrix_sum_test_double("double");

template <typename DataType_>
class DenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(-3)),
                dm3(size, size + 1, DataType_(-1));
            DenseMatrix<DataType_> & sum(Sum<>::value(dm1, dm2));

            TEST_CHECK_EQUAL(sum, dm3);

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(Sum<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSumQuickTest<float> dense_matrix_sum_quick_test_float("float");
DenseMatrixSumQuickTest<double> dense_matrix_sum_quick_test_double("double");

template <typename DataType_>
class SparseMatrixSumTest :
    public BaseTest
{
    public:
        SparseMatrixSumTest(const std::string & type) :
            BaseTest("sparse_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1),
                    sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(2);
                        *k = DataType_(2);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_(-3);
                        *k = DataType_(-3);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = DataType_(-1);
                    }
                }
                SparseMatrix<DataType_> & sum(Sum<>::value(sm1, sm2));

                TEST_CHECK_EQUAL(sum, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Sum<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixSumTest<float> sparse_matrix_sum_test_float("float");
SparseMatrixSumTest<double> sparse_matrix_sum_test_double("double");

template <typename DataType_>
class SparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        SparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1),
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(2);
                    *k = DataType_(2);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_(-3);
                    *k = DataType_(-3);
                }
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                {
                    *k = DataType_(-1);
                }
            }
            SparseMatrix<DataType_> & sum(Sum<>::value(sm1, sm2));

            TEST_CHECK_EQUAL(sum, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Sum<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixSumQuickTest<float> sparse_matrix_sum_quick_test_float("float");
SparseMatrixSumQuickTest<double> sparse_matrix_sum_quick_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixSumTest :
    public BaseTest
{
    public:
        ScalarDenseMatrixSumTest(const std::string & type) :
            BaseTest("scalar_dense_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(5));
                DenseMatrix<DataType_> & prod(Sum<>::value(dm1, DataType_(3)));

                TEST_CHECK_EQUAL(prod, dm2);
            }
        }
};
ScalarDenseMatrixSumTest<float> scalar_dense_matrix_sum_test_float("float");
ScalarDenseMatrixSumTest<double> scalar_dense_matrix_sum_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm(size, size + 1, DataType_(2));
            DenseMatrix<DataType_> prod1(Sum<>::value(dm, DataType_(3)));

            DataType_ vsum(0), ssum(2 * DataType_(size) * DataType_(size + 1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                vsum += *i;
            }
            TEST_CHECK_EQUAL_WITHIN_EPS(vsum, DataType_(5 * size * (size +1)),
                    std::numeric_limits<DataType_>::epsilon());
        }
};
ScalarDenseMatrixSumQuickTest<float> scalar_dense_matrix_sum_quick_test_float("float");
ScalarDenseMatrixSumQuickTest<double> scalar_dense_matrix_sum_quick_test_double("double");

// Test cases for vector operations

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
                DenseVector<DataType_> dv1(size, DataType_(1)), dv2(size, DataType_(2));

                Sum<>::value(dv1, dv2);

                DataType_ v1(Norm<vnt_l_one>::value(dv1));
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, 3 * size, std::numeric_limits<DataType_>::epsilon());

                DenseVector<DataType_> dv3(size + 1), dv4(size);
                TEST_CHECK_THROWS(Sum<>::value(dv3, dv4), VectorSizeDoesNotMatch);
            }
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
            DenseVector<DataType_> dv1(size, DataType_(1)), dv2(size, DataType_(2));

            Sum<>::value(dv1, dv2);
            DataType_ v1(Norm<vnt_l_one>::value(dv1));
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, 3 * size, std::numeric_limits<DataType_>::epsilon());

            DenseVector<DataType_> dv3(6, DataType_(1)), dv4(4, DataType_(1));
            TEST_CHECK_THROWS(Sum<>::value(dv3, dv4), VectorSizeDoesNotMatch);
        }
};
DenseVectorSumQuickTest<float>  dense_vector_sum_quick_test_float("float");
DenseVectorSumQuickTest<double> dense_vector_sum_quick_test_double("double");

template <typename DataType_>
class DenseVectorSparseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSparseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sparse_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv3(size, DataType_(0));
                SparseVector<DataType_> sv2(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(sv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(1);
                        *k = DataType_(1);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_(2);
                        *k = DataType_(2);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = DataType_(3);
                    }
                }

                Sum<>::value(dv1, sv2);

                TEST_CHECK_EQUAL(dv1, dv3);
            }

            DenseVector<DataType_> dv01(2);
            SparseVector<DataType_> sv02(1, 1);
            TEST_CHECK_THROWS(Sum<>::value(dv01, sv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorSumTest<float> dense_vector_sparse_vector_sum_test_float("float");
DenseVectorSparseVectorSumTest<double> dense_vector_sparse_vector_sum_test_double("double");

template <typename DataType_>
class DenseVectorSparseVectorSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSparseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sparse_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> dv1(size, DataType_(0)), dv3(size, DataType_(0));
            SparseVector<DataType_> sv2(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                j(sv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(1);
                    *k = DataType_(1);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_(2);
                    *k = DataType_(2);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    *k = DataType_(3);
                }
            }

            Sum<>::value(dv1, sv2);

            TEST_CHECK_EQUAL(dv1, dv3);

            DenseVector<DataType_> dv01(2);
            SparseVector<DataType_> sv02(1, 1);
            TEST_CHECK_THROWS(Sum<>::value(dv01, sv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorSumQuickTest<float> dense_vector_sparse_vector_sum_quick_test_float("float");
DenseVectorSparseVectorSumQuickTest<double> dense_vector_sparse_vector_sum_quick_test_double("double");

template <typename DataType_>
class SparseVectorSumTest :
    public BaseTest
{
    public:
        SparseVectorSumTest(const std::string & type) :
            BaseTest("sparse_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(1);
                        *k = DataType_(1);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_(2);
                        *k = DataType_(2);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = DataType_(3);
                    }
                }

                Sum<>::value(sv1, sv2);

                TEST_CHECK_EQUAL(sv1, sv3);
            }

            SparseVector<DataType_> sv01(2, 1), sv02(1, 1);
            TEST_CHECK_THROWS(Sum<>::value(sv01, sv02), VectorSizeDoesNotMatch);
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
            QuickTest("sparse_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(1);
                    *k = DataType_(1);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_(2);
                    *k = DataType_(2);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    *k = DataType_(3);
                }
            }

            Sum<>::value(sv1, sv2);

            TEST_CHECK_EQUAL(sv1, sv3);

            SparseVector<DataType_> sv01(21, 1), sv02(20, 1);
            TEST_CHECK_THROWS(Sum<>::value(sv01, sv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorSumQuickTest<float> sparse_vector_sum_quick_test_float("float");
SparseVectorSumQuickTest<double> sparse_vector_sum_quick_test_double("double");
