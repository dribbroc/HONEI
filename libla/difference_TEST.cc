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
#include <libla/difference.hh>
#include <libla/matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;

// Test cases for matrix operations

template <typename DT_>
class BandedMatrixDenseMatrixDifferenceTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixDifferenceTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size, DT_(2));
                BandedMatrix<DT_> bm1(size, dv1);
                DenseMatrix<DT_> dm2(size, size, DT_(1)), dm3(size, size, DT_(-1));

                typename MutableMatrix<DT_>::ElementIterator k(dm3.begin_elements());
                for (typename Matrix<DT_>::ConstElementIterator i(bm1.begin_elements()),
                    i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
                {
                    if (*i != DT_(0))
                    {
                        *k = DT_(1);
                    }
                }
                DenseMatrix<DT_> & difference(Difference<>::value(bm1, dm2));

                TEST_CHECK_EQUAL(difference, dm3);
            }

            BandedMatrix<DT_> bm01(5); 
            DenseMatrix<DT_> dm02(6, 5), dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(bm01, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(bm01, dm02), MatrixColumnsDoNotMatch);
        }
};
BandedMatrixDenseMatrixDifferenceTest<float> banded_matrix_dense_matrix_difference_test_float("float");
BandedMatrixDenseMatrixDifferenceTest<double> banded_matrix_dense_matrix_difference_test_double("double");

template <typename DT_>
class BandedMatrixDenseMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseVector<DT_> dv1(size, DT_(2));
            BandedMatrix<DT_> bm1(size, dv1);
            DenseMatrix<DT_> dm2(size, size, DT_(1)), dm3(size, size, DT_(-1));

            typename MutableMatrix<DT_>::ElementIterator k(dm3.begin_elements());
            for (typename Matrix<DT_>::ConstElementIterator i(bm1.begin_elements()),
                i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
            {
                if (*i != DT_(0))
                {
                    *k = DT_(1);
                }
            }
            DenseMatrix<DT_> & difference(Difference<>::value(bm1, dm2));

            TEST_CHECK_EQUAL(difference, dm3);

            BandedMatrix<DT_> bm01(5); 
            DenseMatrix<DT_> dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(bm01, dm03), MatrixIsNotSquare);
            TEST_CHECK_THROWS(Difference<>::value(bm01, dm02), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseMatrixDifferenceQuickTest<float> banded_matrix_dense_matrix_difference_quick_test_float("float");
BandedMatrixDenseMatrixDifferenceQuickTest<double> banded_matrix_dense_matrix_difference_quick_test_double("double");

template <typename DT_>
class BandedMatrixDifferenceTest :
    public BaseTest
{
    public:
        BandedMatrixDifferenceTest(const std::string & type) :
            BaseTest("banded_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size, DT_(2));
                DenseVector<DT_> dv2(size, DT_(3));
                DenseVector<DT_> dv3(size, DT_(-1));
                BandedMatrix<DT_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);
                BandedMatrix<DT_> & difference(Difference<>::value(bm1, bm2));

                TEST_CHECK_EQUAL(difference, bm3);
            }

            BandedMatrix<DT_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Difference<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixDifferenceTest<float> banded_matrix_difference_test_float("float");
BandedMatrixDifferenceTest<double> banded_matrix_difference_test_double("double");

template <typename DT_>
class BandedMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("banded_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (5);
            DenseVector<DT_> dv1(size, DT_(2)), dv2(size, DT_(3)), dv3(size, DT_(-1));
            BandedMatrix<DT_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);

            BandedMatrix<DT_> & diff(Difference<>::value(bm1, bm2));

            for (typename Matrix<DT_>::ConstElementIterator i(diff.begin_elements()),
                    i_end(diff.end_elements()), j(bm3.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DT_>::epsilon());
            }

            BandedMatrix<DT_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Difference<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixDifferenceQuickTest<float> banded_matrix_difference_quick_test_float("float");
BandedMatrixDifferenceQuickTest<double> banded_matrix_difference_quick_test_double("double");

template <typename DT_>
class BandedMatrixSparseMatrixDifferenceTest :
    public BaseTest
{
    public:
        BandedMatrixSparseMatrixDifferenceTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size, DT_(2));
                BandedMatrix<DT_> bm1(size, dv1);
                SparseMatrix<DT_> sm2(size, size, size / 8 + 1),
                        sm3(size, size, size / 8 + 1);

                for (typename MutableMatrix<DT_>::ElementIterator i(sm2.begin_elements()),
                        i_end(sm2.end_elements()), k(sm3.begin_elements()) ;
                        i != i_end ; ++i, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DT_(1);
                        *k = DT_(-1);
                    }
                }

                typename MutableMatrix<DT_>::ElementIterator k(sm3.begin_elements());
                for (typename Matrix<DT_>::ConstElementIterator i(bm1.begin_elements()),
                        i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
                {
                    if (*i != DT_(0))
                    {
                        *k = DT_(2);
                        if (i.index() % 10 == 0)
                            *k = DT_(1);
                    }
                }
                SparseMatrix<DT_> & difference(Difference<>::value(bm1, sm2));

                TEST_CHECK_EQUAL(difference, sm3);
            }

            BandedMatrix<DT_> bm01(5);
            SparseMatrix<DT_> sm02(6, 5, 1), sm03(6, 6, 1);

            TEST_CHECK_THROWS(Difference<>::value(bm01, sm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(bm01, sm02), MatrixIsNotSquare);
        }
};
BandedMatrixSparseMatrixDifferenceTest<float> banded_matrix_sparse_matrix_difference_test_float("float");
BandedMatrixSparseMatrixDifferenceTest<double> banded_matrix_sparse_matrix_difference_test_double("double");

template <typename DT_>
class BandedMatrixSparseMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DT_> dv1(size, DT_(2));
            BandedMatrix<DT_> bm1(size, dv1);
            SparseMatrix<DT_> sm2(size, size, size / 8 + 1),
                    sm3(size, size, size / 8 + 1);

            for (typename MutableMatrix<DT_>::ElementIterator i(sm2.begin_elements()),
                i_end(sm2.end_elements()), k(sm3.begin_elements()) ;
                i != i_end ; ++i, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DT_(1);
                    *k = DT_(-1);
                }
            }

           typename MutableMatrix<DT_>::ElementIterator k(sm3.begin_elements());
            for (typename Matrix<DT_>::ConstElementIterator i(bm1.begin_elements()),
                i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
            {
                if (*i != DT_(0))
                {
                    *k = DT_(2);
                    if (i.index() % 10 == 0)
                        *k = DT_(1);
                }
            }

            SparseMatrix<DT_> & difference(Difference<>::value(bm1, sm2));
            TEST_CHECK_EQUAL(difference, sm3);

            BandedMatrix<DT_> bm01(5); 
            SparseMatrix<DT_> sm02(6, 5, 1), sm03(6, 6, 1);

            TEST_CHECK_THROWS(Difference<>::value(bm01, sm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(bm01, sm02), MatrixIsNotSquare);
        }
};
BandedMatrixSparseMatrixDifferenceQuickTest<float> banded_matrix_sparse_matrix_difference_quick_test_float("float");
BandedMatrixSparseMatrixDifferenceQuickTest<double> banded_matrix_sparse_matrix_difference_quick_test_double("double");

template <typename DT_>
class DenseMatrixDifferenceTest :
    public BaseTest
{
    public:
        DenseMatrixDifferenceTest(const std::string & type) :
            BaseTest("dense_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size + 1, DT_(2)), dm2(size, size + 1, DT_(3)),
                    dm3(size, size + 1, DT_(-1));
                DenseMatrix<DT_> & difference(Difference<>::value(dm1, dm2));

                TEST_CHECK_EQUAL(difference, dm3);
            }

            DenseMatrix<DT_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixDifferenceTest<float> dense_matrix_difference_test_float("float");
DenseMatrixDifferenceTest<double> dense_matrix_difference_test_double("double");

template <typename DT_>
class DenseMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        DenseMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("dense_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DT_> dm1(size, size + 1, DT_(2)), dm2(size, size + 1, DT_(3)),
                dm3(size, size + 1, DT_(-1));
            DenseMatrix<DT_> & difference(Difference<>::value(dm1, dm2));

            TEST_CHECK_EQUAL(difference, dm3);

            DenseMatrix<DT_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixDifferenceQuickTest<float> dense_matrix_difference_quick_test_float("float");
DenseMatrixDifferenceQuickTest<double> dense_matrix_difference_quick_test_double("double");

template <typename DT_>
class DenseMatrixSparseMatrixDifferenceTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixDifferenceTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size + 1, DT_(0)),
                        dm3(size, size + 1, DT_(0));
                SparseMatrix<DT_> sm2(size, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()),
                    i_end(dm1.end_elements()), j(sm2.begin_elements()), k(dm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DT_((i.index() +1) * 2 / 1.23456789);
                        *k = DT_((i.index() +1) * 2 / 1.23456789);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DT_((i.index() +1) / 1.23456789);
                        *k = DT_((i.index() +1) / -1.23456789);
                    }
                    if (i.index() % 10 == 0 && i.index() % 7 == 0)
                    {
                        *k = DT_((i.index() +1) / 1.23456789);
                    }
                }
                DenseMatrix<DT_> & difference(Difference<>::value(dm1, sm2));

                TEST_CHECK_EQUAL(difference, dm3);
            }

            SparseMatrix<DT_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DT_> dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixDifferenceTest<float> dense_matrix_sparse_matrix_difference_test_float("float");
DenseMatrixSparseMatrixDifferenceTest<double> dense_matrix_sparse_matrix_difference_test_double("double");

template <typename DT_>
class DenseMatrixSparseMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (10);
            DenseMatrix<DT_> dm1(size, size + 1, DT_(0)),
                    dm3(size, size + 1, DT_(0));
            SparseMatrix<DT_> sm2(size, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()),
                i_end(dm1.end_elements()), j(sm2.begin_elements()), k(dm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DT_((i.index() +1) * 2 / 1.23456789);
                    *k = DT_((i.index() +1) * 2 / 1.23456789);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DT_((i.index() +1) / 1.23456789);
                    *k = DT_((i.index() +1) / -1.23456789);
                }
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                {
                    *k = DT_((i.index() +1) / 1.23456789);
                }
            }
            DenseMatrix<DT_> & difference(Difference<>::value(dm1, sm2));

            TEST_CHECK_EQUAL(difference, dm3);

            SparseMatrix<DT_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DT_> dm03(5, 6);

            TEST_CHECK_THROWS(Difference<>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixDifferenceQuickTest<float> dense_matrix_sparse_matrix_difference_quick_test_float("float");
DenseMatrixSparseMatrixDifferenceQuickTest<double> dense_matrix_sparse_matrix_difference_quick_test_double("double");

template <typename DT_>
class SparseMatrixDifferenceTest :
    public BaseTest
{
    public:
        SparseMatrixDifferenceTest(const std::string & type) :
            BaseTest("sparse_matrix_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size + 1, size / 8 + 1),
                    sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DT_((i.index() +1) * 2 / 1.23456789);
                        *k = DT_((i.index() +1) * 2 / 1.23456789);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DT_((i.index() +1) / 1.23456789);
                        *k = DT_((i.index() +1) / -1.23456789);
                    }
                    if (i.index() % 10 == 0 && i.index() % 7 == 0)
                    {
                        *k = DT_((i.index() +1) / 1.23456789);
                    }
                }
                SparseMatrix<DT_> & difference(Difference<>::value(sm1, sm2));

                TEST_CHECK_EQUAL(difference, sm3);
            }

            SparseMatrix<DT_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Difference<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDifferenceTest<float> sparse_matrix_difference_test_float("float");
SparseMatrixDifferenceTest<double> sparse_matrix_difference_test_double("double");

template <typename DT_>
class SparseMatrixDifferenceQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDifferenceQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            SparseMatrix<DT_> sm1(size, size + 1, size / 8 + 1),
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DT_((i.index() +1) * 2 / 1.23456789);
                    *k = DT_((i.index() +1) * 2 / 1.23456789);
                }
                if (i.index() % 7 == 0)
                {
                    *j = DT_((i.index() +1) / 1.23456789);
                    *k = DT_((i.index() +1) / -1.23456789);
                }
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                {
                    *k = DT_((i.index() +1) / 1.23456789);
                }
            }
            SparseMatrix<DT_> & difference(Difference<>::value(sm1, sm2));

            TEST_CHECK_EQUAL(difference, sm3);

            SparseMatrix<DT_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Difference<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Difference<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDifferenceQuickTest<float> sparse_matrix_difference_quick_test_float("float");
SparseMatrixDifferenceQuickTest<double> sparse_matrix_difference_quick_test_double("double");

// Test cases for vector operations

template <typename DT_>
class DenseVectorDifferenceTest :
    public BaseTest
{
    public:
        DenseVectorDifferenceTest(const std::string & type) :
            BaseTest("dense_vector_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size), dv2(size), dv3(size, DT_(0));

                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
                }
                for (typename Vector<DT_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
                }
                DenseVector<DT_> difference1(Difference<>::value(dv1, dv2));
                TEST_CHECK_EQUAL(difference1, dv3);
            }

            DenseVector<DT_> dv00(1, DT_(1));
            DenseVector<DT_> dv01(5, DT_(1));
            TEST_CHECK_THROWS(Difference<DT_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorDifferenceTest<float> dense_vector_difference_test_float("float");
DenseVectorDifferenceTest<double> dense_vector_difference_test_double("double");

template <typename DT_>
class DenseVectorDifferenceQuickTest :
    public QuickTest
{
    public:
        DenseVectorDifferenceQuickTest(const std::string & type) :
            QuickTest("dense_vector_difference_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DT_> dv1(size), dv2(size), dv3(size, DT_(0));

            for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
            }
            for (typename Vector<DT_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
            }
            DenseVector<DT_> difference1(Difference<>::value(dv1, dv2));
            TEST_CHECK_EQUAL(difference1, dv3);

            DenseVector<DT_> dv00(1, DT_(1));
            DenseVector<DT_> dv01(5, DT_(1));

            TEST_CHECK_THROWS(Difference<DT_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorDifferenceQuickTest<float>  dense_vector_difference_quick_test_float("float");
DenseVectorDifferenceQuickTest<double> dense_vector_difference_quick_test_double("double");

template <typename DT_>
class DenseVectorSparseVectorDifferenceTest :
    public BaseTest
{
    public:
        DenseVectorSparseVectorDifferenceTest(const std::string & type) :
            BaseTest("dense_vector_sparse_vector_difference_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size, DT_(0)), dv3(size, DT_(0));
                SparseVector<DT_> sv2(size, size / 7 + 1);
                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(sv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = static_cast<DT_>(((i.index() +1) * 2) / 1.23456789);
                        *k = static_cast<DT_>(((i.index() +1) * 2) / 1.23456789);
                    }
                    if (i.index() % 7 == 0) 
                    {
                        *j = static_cast<DT_>((i.index() +1) / 1.23456789);
                        *k = static_cast<DT_>((i.index() +1) / -1.23456789);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = static_cast<DT_>((i.index() +1) / 1.23456789);
                    }
                }

                DenseVector<DT_> difference1(Difference<>::value(dv1, sv2));
                TEST_CHECK_EQUAL(difference1, dv3);
            }

            DenseVector<DT_> dv00(1);
            SparseVector<DT_> sv01(5, 1);
            TEST_CHECK_THROWS(Difference<DT_>::value(dv00, sv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorDifferenceTest<float> dense_vector_sparse_vector_difference_test_float("float");
DenseVectorSparseVectorDifferenceTest<double> dense_vector_sparse_vector_difference_test_double("double");

template <typename DT_>
class DenseVectorSparseVectorDifferenceQuickTest :
    public QuickTest
{
    public:
        DenseVectorSparseVectorDifferenceQuickTest(const std::string & type) :
            QuickTest("dense_vector_sparse_vector_difference_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DT_> dv1(size, DT_(0)), dv3(size, DT_(0));
            SparseVector<DT_> sv2(size, size / 7 + 1);
            for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(sv2.begin_elements()), k(dv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = static_cast<DT_>(((i.index() +1) * 2) / 1.23456789);
                    *k = static_cast<DT_>(((i.index() +1) * 2) / 1.23456789);
                }
                if (i.index() % 7 == 0) 
                {
                    *j = static_cast<DT_>((i.index() +1) / 1.23456789);
                    *k = static_cast<DT_>((i.index() +1) / -1.23456789);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    *k = static_cast<DT_>((i.index() +1) / 1.23456789);
                }
            }

            DenseVector<DT_> difference1(Difference<>::value(dv1, sv2));
            TEST_CHECK_EQUAL(difference1, dv3);

            DenseVector<DT_> dv00(1);
            SparseVector<DT_> sv01(5, 1);
            TEST_CHECK_THROWS(Difference<DT_>::value(dv00, sv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorDifferenceQuickTest<float> dense_vector_sparse_vector_difference_quick_test_float("float");
DenseVectorSparseVectorDifferenceQuickTest<double> dense_vector_sparse_vector_difference_quick_test_double("double");

template <typename DT_>
class SparseVectorDifferenceTest :
    public BaseTest
{
    public:
        SparseVectorDifferenceTest(const std::string & type) :
            BaseTest("sparse_vector_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DT_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if ((i.index() % 7 == 0) && (i.index() % 10 == 0))
                    {
                        *i = DT_((i.index() + 1) / 1.23456789);
                        *j = DT_(((i.index() +1) * 2) / 1.23456789);
                        *k = -*i;
                    }
                    else if (i.index() % 7 == 0) 
                    {
                        *i = DT_((i.index() +1) / 1.23456789);
                        *k = *i;
                    }
                    else if (i.index() % 10 == 0) 
                    {
                        *j = DT_(((i.index() +1) * 2) / 1.23456789);
                        *k = -*j;
                    }
                }

                SparseVector<DT_> difference1(Difference<>::value(sv1, sv2));
                TEST_CHECK_EQUAL(difference1, sv3);
            }

            SparseVector<DT_> sv00(1, 1), sv01(5, 1);
            TEST_CHECK_THROWS(Difference<DT_>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorDifferenceTest<float> sparse_vector_difference_test_float("float");
SparseVectorDifferenceTest<double> sparse_vector_difference_test_double("double");

template <typename DT_>
class SparseVectorDifferenceQuickTest :
    public QuickTest
{
    public:
        SparseVectorDifferenceQuickTest(const std::string & type) :
            QuickTest("sparse_vector_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(25);
            SparseVector<DT_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
            for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if ((i.index() % 7 == 0) && (i.index() % 10 == 0))
                {
                    *i = DT_((i.index() + 1) / 1.23456789);
                    *j = DT_(((i.index() +1) * 2) / 1.23456789);
                    *k = -*i;
                }
                else if (i.index() % 7 == 0) 
                {
                    *i = DT_((i.index() +1) / 1.23456789);
                    *k = *i;
                }
                else if (i.index() % 10 == 0) 
                {
                    *j = DT_(((i.index() +1) * 2) / 1.23456789);
                    *k = -*j;
                }
            }

            SparseVector<DT_> difference1(Difference<>::value(sv1, sv2));
            TEST_CHECK_EQUAL(difference1, sv3);

            SparseVector<DT_> sv00(1, 1), sv01(5, 1);
            TEST_CHECK_THROWS(Difference<DT_>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorDifferenceQuickTest<float> sparse_vector_difference_quick_test_float("float");
SparseVectorDifferenceQuickTest<double> sparse_vector_difference_quick_test_double("double");

