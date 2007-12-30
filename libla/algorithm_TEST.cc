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
#include <libla/algorithm.hh>
#include <unittest/unittest.hh>

using namespace honei;
using namespace tests;


class BandedMatrixConvertTest :
    public BaseTest
{
    public:
        BandedMatrixConvertTest() :
            BaseTest("banded_matrix_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 6) ; size <<= 1)
            {
                BandedMatrix<float> bmf1(size), bmf2(size), bmfr(size);
                BandedMatrix<double> bmd1(size), bmd2(size), bmdr(size);
                MutableMatrix<float>::ElementIterator fr(bmfr.begin_elements());
                MutableMatrix<double>::ElementIterator dr(bmdr.begin_elements());
                MutableMatrix<float>::ElementIterator f1(bmf1.begin_elements());
                MutableMatrix<double>::ElementIterator d1(bmd1.begin_elements());
                for (MutableMatrix<float>::ElementIterator i_end(bmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    if (fr.index() % 50 == 0)
                    {
                        *fr = float(fr.index()) / 1.1234;
                        *dr = double(fr.index()) / 1.1234;
                        *f1 = float(fr.index()) / 1.1234;
                        *d1 = double(fr.index()) / 1.1234;
                    }
                }
                convert(bmf2, bmd1);
                convert(bmd2, bmf1);
                TEST_CHECK_EQUAL(bmf2, bmfr);
                for (Matrix<double>::ConstElementIterator i(bmd2.begin_elements()), r(bmdr.begin_elements()), i_end(bmd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
            }
            BandedMatrix<float> bmf01(3), bmf02(3);
            BandedMatrix<double> bmd01(2), bmd02(4);
            TEST_CHECK_THROWS(convert(bmf01, bmd01), MatrixSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(bmd02, bmf02), MatrixSizeDoesNotMatch);
        }
} banded_matrix_convert_test;

class DenseMatrixConvertTest :
    public BaseTest
{
    public:
        DenseMatrixConvertTest() :
            BaseTest("dense_matrix_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 6) ; size <<= 1)
            {
                DenseMatrix<float> dmf1(size, size + 1), dmf2(size, size + 1), dmfr(size, size + 1);
                DenseMatrix<double> dmd1(size, size + 1), dmd2(size, size + 1), dmdr(size, size + 1);
                MutableMatrix<float>::ElementIterator fr(dmfr.begin_elements());
                MutableMatrix<double>::ElementIterator dr(dmdr.begin_elements());
                MutableMatrix<float>::ElementIterator f1(dmf1.begin_elements());
                MutableMatrix<double>::ElementIterator d1(dmd1.begin_elements());
                for (MutableMatrix<float>::ElementIterator i_end(dmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
                convert(dmf2, dmd1);
                convert(dmd2, dmf1);
                TEST_CHECK_EQUAL(dmf2, dmfr);
                for (Matrix<double>::ConstElementIterator i(dmd2.begin_elements()), r(dmdr.begin_elements()), i_end(dmd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
            }
            DenseMatrix<float> dmf01(3, 2), dmf02(3, 3);
            DenseMatrix<double> dmd01(2, 2), dmd02(3, 4);
            TEST_CHECK_THROWS(convert(dmf01, dmd01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(convert(dmd02, dmf02), MatrixColumnsDoNotMatch);
        }
} dense_matrix_convert_test;

class DenseMatrixConvertQuickTest :
    public QuickTest
{
    public:
        DenseMatrixConvertQuickTest() :
            QuickTest("dense_matrix_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(47);
            DenseMatrix<float> dmf1(size, size + 1), dmf2(size, size + 1), dmfr(size, size + 1);
            DenseMatrix<double> dmd1(size, size + 1), dmd2(size, size + 1), dmdr(size, size + 1);
            MutableMatrix<float>::ElementIterator fr(dmfr.begin_elements());
            MutableMatrix<double>::ElementIterator dr(dmdr.begin_elements());
            MutableMatrix<float>::ElementIterator f1(dmf1.begin_elements());
            MutableMatrix<double>::ElementIterator d1(dmd1.begin_elements());
            for (MutableMatrix<float>::ElementIterator i_end(dmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                *fr = float(fr.index()) / 1.1234;
                *dr = double(fr.index()) / 1.1234;
                *f1 = float(fr.index()) / 1.1234;
                *d1 = double(fr.index()) / 1.1234;
            }
            convert(dmf2, dmd1);
            convert(dmd2, dmf1);
            TEST_CHECK_EQUAL(dmf2, dmfr);
            for (Matrix<double>::ConstElementIterator i(dmd2.begin_elements()), r(dmdr.begin_elements()), i_end(dmd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            DenseMatrix<float> dmf01(3, 2), dmf02(3, 3);
            DenseMatrix<double> dmd01(2, 2), dmd02(3, 4);
            TEST_CHECK_THROWS(convert(dmf01, dmd01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(convert(dmd02, dmf02), MatrixColumnsDoNotMatch);
        }
} dense_matrix_convert_quick_test;

class DenseVectorConvertTest :
    public BaseTest
{
    public:
        DenseVectorConvertTest() :
            BaseTest("dense_vector_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<float> dvf1(size), dvf2(size), dvfr(size);
                DenseVector<double> dvd1(size), dvd2(size), dvdr(size);
                Vector<float>::ElementIterator fr(dvfr.begin_elements());
                Vector<double>::ElementIterator dr(dvdr.begin_elements());
                Vector<float>::ElementIterator f1(dvf1.begin_elements());
                Vector<double>::ElementIterator d1(dvd1.begin_elements());
                for (Vector<float>::ElementIterator i_end(dvf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
                convert(dvf2, dvd1);
                convert(dvd2, dvf1);
                TEST_CHECK_EQUAL(dvf2, dvfr);
                for (Vector<double>::ConstElementIterator i(dvd2.begin_elements()), r(dvdr.begin_elements()), i_end(dvd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
            }
            DenseVector<float> dvf01(3), dvf02(3);
            DenseVector<double> dvd01(2), dvd02(4);
            TEST_CHECK_THROWS(convert(dvf01, dvd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(dvd02, dvf02), VectorSizeDoesNotMatch);
        }
} dense_vector_convert_test;

class DenseVectorConvertQuickTest :
    public QuickTest
{
    public:
        DenseVectorConvertQuickTest() :
            QuickTest("dense_vector_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(4711);
            DenseVector<float> dvf1(size), dvf2(size), dvfr(size);
            DenseVector<double> dvd1(size), dvd2(size), dvdr(size);
            Vector<float>::ElementIterator fr(dvfr.begin_elements());
            Vector<double>::ElementIterator dr(dvdr.begin_elements());
            Vector<float>::ElementIterator f1(dvf1.begin_elements());
            Vector<double>::ElementIterator d1(dvd1.begin_elements());
            for (Vector<float>::ElementIterator i_end(dvf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                *fr = float(fr.index()) / 1.1234;
                *dr = double(fr.index()) / 1.1234;
                *f1 = float(fr.index()) / 1.1234;
                *d1 = double(fr.index()) / 1.1234;
            }
            convert(dvf2, dvd1);
            convert(dvd2, dvf1);
            TEST_CHECK_EQUAL(dvf2, dvfr);
            for (Vector<double>::ConstElementIterator i(dvd2.begin_elements()), r(dvdr.begin_elements()), i_end(dvd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            DenseVector<float> dvf01(3), dvf02(3);
            DenseVector<double> dvd01(2), dvd02(4);
            TEST_CHECK_THROWS(convert(dvf01, dvd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(dvd02, dvf02), VectorSizeDoesNotMatch);
        }
} dense_vector_convert_quick_test;

class DenseVectorRangeConvertTest :
    public BaseTest
{
    public:
        DenseVectorRangeConvertTest() :
            BaseTest("dense_vector_range_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                for (int a(0) ; a < 4 ; a++)
                {
                    for (int b(0) ; b < 4 ; b++)
                    {
                        DenseVector<float> dvf1(size * 2), dvf2(size * 2), dvfr(size * 2);
                        DenseVector<double> dvd1(size * 2), dvd2(size * 2), dvdr(size * 2);
                        DenseVectorRange<float> dvf1r(dvf1, size, a), dvf2r(dvf2, size, a), dvfrr(dvfr, size, a);
                        DenseVectorRange<double> dvd1r(dvd1, size, b), dvd2r(dvd2, size, b), dvdrr(dvdr, size, b);
                        Vector<float>::ElementIterator fr(dvfrr.begin_elements());
                        Vector<double>::ElementIterator dr(dvdrr.begin_elements());
                        Vector<float>::ElementIterator f1(dvf1r.begin_elements());
                        Vector<double>::ElementIterator d1(dvd1r.begin_elements());
                        for (Vector<float>::ElementIterator i_end(dvf1r.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                        {
                            *fr = float(fr.index()) / 1.1234;
                            *dr = double(fr.index()) / 1.1234;
                            *f1 = float(fr.index()) / 1.1234;
                            *d1 = double(fr.index()) / 1.1234;
                        }
                        convert(dvf2r, dvd1r);
                        convert(dvd2r, dvf1r);
                        TEST_CHECK_EQUAL(dvf2r, dvfrr);
                        for (Vector<double>::ConstElementIterator i(dvd2r.begin_elements()), r(dvdrr.begin_elements()), i_end(dvd2r.end_elements())
                                ; i != i_end ; ++i, ++r)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                        }
                    }
                }
            }
            DenseVector<float> dvf01(3), dvf02(3);
            DenseVector<double> dvd01(2), dvd02(4);
            TEST_CHECK_THROWS(convert(dvf01, dvd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(dvd02, dvf02), VectorSizeDoesNotMatch);
        }
} dense_vector_range_convert_test;

class DenseVectorRangeConvertQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeConvertQuickTest() :
            QuickTest("dense_vector_range_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(1711);
            for (int a(0) ; a < 4 ; a++)
            {
                for (int b(0) ; b < 4 ; b++)
                {
                    DenseVector<float> dvf1(size * 2), dvf2(size * 2), dvfr(size * 2);
                    DenseVector<double> dvd1(size * 2), dvd2(size * 2), dvdr(size * 2);
                    DenseVectorRange<float> dvf1r(dvf1, size, a), dvf2r(dvf2, size, a), dvfrr(dvfr, size, a);
                    DenseVectorRange<double> dvd1r(dvd1, size, b), dvd2r(dvd2, size, b), dvdrr(dvdr, size, b);
                    Vector<float>::ElementIterator fr(dvfrr.begin_elements());
                    Vector<double>::ElementIterator dr(dvdrr.begin_elements());
                    Vector<float>::ElementIterator f1(dvf1r.begin_elements());
                    Vector<double>::ElementIterator d1(dvd1r.begin_elements());
                    for (Vector<float>::ElementIterator i_end(dvf1r.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                    {
                        *fr = float(fr.index()) / 1.1234;
                        *dr = double(fr.index()) / 1.1234;
                        *f1 = float(fr.index()) / 1.1234;
                        *d1 = double(fr.index()) / 1.1234;
                    }
                    convert(dvf2r, dvd1r);
                    convert(dvd2r, dvf1r);
                    TEST_CHECK_EQUAL(dvf2r, dvfrr);
                    for (Vector<double>::ConstElementIterator i(dvd2r.begin_elements()), r(dvdrr.begin_elements()), i_end(dvd2r.end_elements())
                            ; i != i_end ; ++i, ++r)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                    }
                }
            }
            DenseVector<float> dvf01(3), dvf02(3);
            DenseVector<double> dvd01(2), dvd02(4);
            TEST_CHECK_THROWS(convert(dvf01, dvd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(dvd02, dvf02), VectorSizeDoesNotMatch);
        }
} dense_vector_range_convert_quick_test;

class SparseMatrixConvertTest :
    public BaseTest
{
    public:
        SparseMatrixConvertTest() :
            BaseTest("sparse_matrix_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 6) ; size <<= 1)
            {
                SparseMatrix<float> smf1(size, size + 1), smf2(size, size + 1), smfr(size, size + 1);
                SparseMatrix<double> smd1(size, size + 1), smd2(size, size + 1), smdr(size, size + 1);
                MutableMatrix<float>::ElementIterator fr(smfr.begin_elements());
                MutableMatrix<double>::ElementIterator dr(smdr.begin_elements());
                MutableMatrix<float>::ElementIterator f1(smf1.begin_elements());
                MutableMatrix<double>::ElementIterator d1(smd1.begin_elements());
                for (MutableMatrix<float>::ElementIterator i_end(smf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    if (fr.index() % 15 == 0)
                    {
                        *fr = float(fr.index()) / 1.1234;
                        *dr = double(fr.index()) / 1.1234;
                        *f1 = float(fr.index()) / 1.1234;
                        *d1 = double(fr.index()) / 1.1234;
                    }
                }
                convert(smf2, smd1);
                convert(smd2, smf1);
                TEST_CHECK_EQUAL(smf2, smfr);
                for (Matrix<double>::ConstElementIterator i(smd2.begin_elements()), r(smdr.begin_elements()), i_end(smd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
            }
            SparseMatrix<float> smf01(3, 2, 1), smf02(3, 3, 1);
            SparseMatrix<double> smd01(2, 2, 1), smd02(3, 4, 1);
            TEST_CHECK_THROWS(convert(smf01, smd01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(convert(smd02, smf02), MatrixColumnsDoNotMatch);
        }
} sparse_matrix_convert_test;

class SparseMatrixConvertQuickTest :
    public QuickTest
{
    public:
        SparseMatrixConvertQuickTest() :
            QuickTest("sparse_matrix_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(47);
            SparseMatrix<float> smf1(size, size + 1), smf2(size, size + 1), smfr(size, size + 1);
            SparseMatrix<double> smd1(size, size + 1), smd2(size, size + 1), smdr(size, size + 1);
            MutableMatrix<float>::ElementIterator fr(smfr.begin_elements());
            MutableMatrix<double>::ElementIterator dr(smdr.begin_elements());
            MutableMatrix<float>::ElementIterator f1(smf1.begin_elements());
            MutableMatrix<double>::ElementIterator d1(smd1.begin_elements());
            for (MutableMatrix<float>::ElementIterator i_end(smf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                if (fr.index() % 15 == 0)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
            }
            convert(smf2, smd1);
            convert(smd2, smf1);
            TEST_CHECK_EQUAL(smf2, smfr);
            for (Matrix<double>::ConstElementIterator i(smd2.begin_elements()), r(smdr.begin_elements()), i_end(smd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            SparseMatrix<float> smf01(3, 2, 1), smf02(3, 3, 1);
            SparseMatrix<double> smd01(2, 2, 1), smd02(3, 4, 1);
            TEST_CHECK_THROWS(convert(smf01, smd01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(convert(smd02, smf02), MatrixColumnsDoNotMatch);
        }
} sparse_matrix_convert_quick_test;

class SparseVectorConvertTest :
    public BaseTest
{
    public:
        SparseVectorConvertTest() :
            BaseTest("sparse_vector_convert_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                SparseVector<float> svf1(size, size / 2), svf2(size, size / 2), svfr(size, size / 2);
                SparseVector<double> svd1(size, size / 2), svd2(size, size / 2), svdr(size, size / 2);
                Vector<float>::ElementIterator fr(svfr.begin_elements());
                Vector<double>::ElementIterator dr(svdr.begin_elements());
                Vector<float>::ElementIterator f1(svf1.begin_elements());
                Vector<double>::ElementIterator d1(svd1.begin_elements());
                for (Vector<float>::ElementIterator i_end(svf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    if (fr.index() % 15 == 0)
                    {
                        *fr = float(fr.index()) / 1.1234;
                        *dr = double(fr.index()) / 1.1234;
                        *f1 = float(fr.index()) / 1.1234;
                        *d1 = double(fr.index()) / 1.1234;
                    }
                }
                convert(svf2, svd1);
                convert(svd2, svf1);
                TEST_CHECK_EQUAL(svf2, svfr);
                for (Vector<double>::ConstElementIterator i(svd2.begin_elements()), r(svdr.begin_elements()), i_end(svd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
            }
            SparseVector<float> svf01(3, 1), svf02(3, 1);
            SparseVector<double> svd01(2, 1), svd02(4, 1);
            TEST_CHECK_THROWS(convert(svf01, svd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(svd02, svf02), VectorSizeDoesNotMatch);
        }
} sparse_vector_convert_test;

class SparseVectorConvertQuickTest :
    public QuickTest
{
    public:
        SparseVectorConvertQuickTest() :
            QuickTest("sparse_vector_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(4711);
            SparseVector<float> svf1(size, size / 2), svf2(size, size / 2), svfr(size, size / 2);
            SparseVector<double> svd1(size, size / 2), svd2(size, size / 2), svdr(size, size / 2);
            Vector<float>::ElementIterator fr(svfr.begin_elements());
            Vector<double>::ElementIterator dr(svdr.begin_elements());
            Vector<float>::ElementIterator f1(svf1.begin_elements());
            Vector<double>::ElementIterator d1(svd1.begin_elements());
            for (Vector<float>::ElementIterator i_end(svf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                if (fr.index() % 15 == 0)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
            }
            convert(svf2, svd1);
            convert(svd2, svf1);
            TEST_CHECK_EQUAL(svf2, svfr);
            for (Vector<double>::ConstElementIterator i(svd2.begin_elements()), r(svdr.begin_elements()), i_end(svd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            SparseVector<float> svf01(3, 1), svf02(3, 1);
            SparseVector<double> svd01(2, 1), svd02(4, 1);
            TEST_CHECK_THROWS(convert(svf01, svd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(svd02, svf02), VectorSizeDoesNotMatch);
        }
} sparse_vector_convert_quick_test;
