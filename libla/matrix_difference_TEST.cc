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
#include <libla/matrix_difference.hh>
#include <libla/matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
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
                DenseVector<DataType_> * dv1(new DenseVector<DataType_>(size, DataType_(2)));
                DenseVector<DataType_> * dv2(new DenseVector<DataType_>(size, DataType_(3)));
                DenseVector<DataType_> * dv3(new DenseVector<DataType_>(size, DataType_(-1)));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);
                BandedMatrix<DataType_> & difference(MatrixDifference<>::value(bm1, bm2));

                TEST_CHECK_EQUAL(difference, bm3);
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixDifference<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixDifferenceTest<float> banded_matrix_difference_test_float("float");
BandedMatrixDifferenceTest<double> banded_matrix_difference_test_double("double");

template <typename DataType_>
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
            DenseVector<DataType_> * dv1(new DenseVector<DataType_>(size, DataType_(2))),
                * dv2(new DenseVector<DataType_>(size, DataType_(3))),
                * dv3(new DenseVector<DataType_>(size, DataType_(-1)));
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);

            BandedMatrix<DataType_> & diff(MatrixDifference<>::value(bm1, bm2));

            for (typename Matrix<DataType_>::ConstElementIterator i(diff.begin_elements()),
                    i_end(diff.end_elements()), j(bm3.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixDifference<>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixDifferenceQuickTest<float> banded_matrix_difference_quick_test_float("float");
BandedMatrixDifferenceQuickTest<double> banded_matrix_difference_quick_test_double("double");

template <typename DataType_>
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
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(3)),
                    dm3(size, size + 1, DataType_(-1));
                DenseMatrix<DataType_> & difference(MatrixDifference<>::value(dm1, dm2));

                TEST_CHECK_EQUAL(difference, dm3);
            }

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(MatrixDifference<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixDifference<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixDifferenceTest<float> dense_matrix_difference_test_float("float");
DenseMatrixDifferenceTest<double> dense_matrix_difference_test_double("double");

template <typename DataType_>
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
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(3)),
                dm3(size, size + 1, DataType_(-1));
            DenseMatrix<DataType_> & difference(MatrixDifference<>::value(dm1, dm2));

            TEST_CHECK_EQUAL(difference, dm3);

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(MatrixDifference<>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixDifference<>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixDifferenceQuickTest<float> dense_matrix_difference_quick_test_float("float");
DenseMatrixDifferenceQuickTest<double> dense_matrix_difference_quick_test_double("double");

template <typename DataType_>
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
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1),
                    sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_((i.index() +1) * 2 / 1.23456789);
                        *k = DataType_((i.index() +1) * 2 / 1.23456789);                        
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = DataType_((i.index() +1) / 1.23456789);
                        *k = DataType_((i.index() +1) / -1.23456789);
                    }
                    if (i.index() % 10 == 0 && i.index() % 7 == 0)
                    {
                        *k = DataType_((i.index() +1) / 1.23456789);
                    }                                        
                }
                SparseMatrix<DataType_> & difference(MatrixDifference<>::value(sm1, sm2));

                TEST_CHECK_EQUAL(difference, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixDifference<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixDifference<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDifferenceTest<float> sparse_matrix_difference_test_float("float");
SparseMatrixDifferenceTest<double> sparse_matrix_difference_test_double("double");

template <typename DataType_>
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
            unsigned long size (21);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1),
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ;
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_((i.index() +1) * 2 / 1.23456789);
                    *k = DataType_((i.index() +1) * 2 / 1.23456789);                        
                }
                if (i.index() % 7 == 0)
                {
                    *j = DataType_((i.index() +1) / 1.23456789);
                    *k = DataType_((i.index() +1) / -1.23456789);
                }
                if (i.index() % 10 == 0 && i.index() % 7 == 0)
                {
                    *k = DataType_((i.index() +1) / 1.23456789);
                }                                        
            }
            SparseMatrix<DataType_> & difference(MatrixDifference<>::value(sm1, sm2));

            TEST_CHECK_EQUAL(difference, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixDifference<>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixDifference<>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixDifferenceQuickTest<float> sparse_matrix_difference_quick_test_float("float");
SparseMatrixDifferenceQuickTest<double> sparse_matrix_difference_quick_test_double("double");
