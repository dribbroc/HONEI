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

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/algorithm.hh>
#include <honei/util/unittest.hh>

#include <iostream>

using namespace honei;
using namespace tests;

template <typename DT_>
class IteratorToDenseVectorCopyTest :
    public QuickTest
{
    public:
        IteratorToDenseVectorCopyTest(const std::string & type) :
            QuickTest("iterator_to_dense_vector_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DT_> dm(17, 9);
            for (typename DenseMatrix<DT_>::ElementIterator i(dm.begin_elements()), i_end(dm.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DT_(i.row() + i.column() * 0.1f);
            }

            DenseVectorSlice<DT_> dvs(dm.column(0));

            DenseVector<DT_> dv(16);
            typename DenseVectorSlice<DT_>::ElementIterator begin(++dvs.begin_elements()), end(dvs.end_elements());
            honei::copy(begin, end, dv);

            for (typename DenseVector<DT_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                    i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, (i.index() + 1), std::numeric_limits<DT_>::epsilon());
            }
        }
};
IteratorToDenseVectorCopyTest<float> iterator_to_dense_vector_copy_test_float("float");
IteratorToDenseVectorCopyTest<double> iterator_to_dense_vector_copy_test_double("double");

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
                BandedMatrix<float>::ElementIterator fr(bmfr.begin_elements());
                BandedMatrix<double>::ElementIterator dr(bmdr.begin_elements());
                BandedMatrix<float>::ElementIterator f1(bmf1.begin_elements());
                BandedMatrix<double>::ElementIterator d1(bmd1.begin_elements());
                for (BandedMatrix<float>::ElementIterator i_end(bmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
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
                for (BandedMatrix<double>::ConstElementIterator i(bmd2.begin_elements()), r(bmdr.begin_elements()), i_end(bmd2.end_elements())
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

class BandedMatrixConvertQuickTest :
    public QuickTest
{
    public:
        BandedMatrixConvertQuickTest() :
            QuickTest("banded_matrix_convert_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned long size(47);
            BandedMatrix<float> bmf1(size), bmf2(size), bmfr(size);
            BandedMatrix<double> bmd1(size), bmd2(size), bmdr(size);
            BandedMatrix<float>::ElementIterator fr(bmfr.begin_elements());
            BandedMatrix<double>::ElementIterator dr(bmdr.begin_elements());
            BandedMatrix<float>::ElementIterator f1(bmf1.begin_elements());
            BandedMatrix<double>::ElementIterator d1(bmd1.begin_elements());
            for (BandedMatrix<float>::ElementIterator i_end(bmf1.end_elements()) ; f1 < i_end ; ++fr, ++dr, ++f1, ++d1)
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
            for (BandedMatrix<double>::ConstElementIterator i(bmd2.begin_elements()), r(bmdr.begin_elements()), i_end(bmd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            BandedMatrix<float> bmf01(3), bmf02(3);
            BandedMatrix<double> bmd01(2), bmd02(4);
            TEST_CHECK_THROWS(convert(bmf01, bmd01), MatrixSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(bmd02, bmf02), MatrixSizeDoesNotMatch);
        }
} banded_matrix_convert_quick_test;

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
                DenseMatrix<float>::ElementIterator fr(dmfr.begin_elements());
                DenseMatrix<double>::ElementIterator dr(dmdr.begin_elements());
                DenseMatrix<float>::ElementIterator f1(dmf1.begin_elements());
                DenseMatrix<double>::ElementIterator d1(dmd1.begin_elements());
                for (DenseMatrix<float>::ElementIterator i_end(dmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
                convert(dmf2, dmd1);
                convert(dmd2, dmf1);
                TEST_CHECK_EQUAL(dmf2, dmfr);
                for (DenseMatrix<double>::ConstElementIterator i(dmd2.begin_elements()), r(dmdr.begin_elements()), i_end(dmd2.end_elements())
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
            DenseMatrix<float>::ElementIterator fr(dmfr.begin_elements());
            DenseMatrix<double>::ElementIterator dr(dmdr.begin_elements());
            DenseMatrix<float>::ElementIterator f1(dmf1.begin_elements());
            DenseMatrix<double>::ElementIterator d1(dmd1.begin_elements());
            for (DenseMatrix<float>::ElementIterator i_end(dmf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                *fr = float(fr.index()) / 1.1234;
                *dr = double(fr.index()) / 1.1234;
                *f1 = float(fr.index()) / 1.1234;
                *d1 = double(fr.index()) / 1.1234;
            }
            convert(dmf2, dmd1);
            convert(dmd2, dmf1);
            TEST_CHECK_EQUAL(dmf2, dmfr);
            for (DenseMatrix<double>::ConstElementIterator i(dmd2.begin_elements()), r(dmdr.begin_elements()), i_end(dmd2.end_elements())
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

template <typename Tag_>
class DenseVectorConvertTest :
    public BaseTest
{
    public:
        DenseVectorConvertTest() :
            BaseTest("dense_vector_convert_test")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<float> dvf1(size), dvf2(size), dvfr(size), dvf3(size), dvff(size);
                DenseVector<double> dvd1(size), dvd2(size), dvdr(size), dvd3(size), dvdf(size);
                DenseVector<float>::ElementIterator fr(dvfr.begin_elements());
                DenseVector<double>::ElementIterator dr(dvdr.begin_elements());
                DenseVector<float>::ElementIterator f1(dvf1.begin_elements());
                DenseVector<double>::ElementIterator d1(dvd1.begin_elements());
                for (DenseVector<float>::ElementIterator i_end(dvf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                {
                    *fr = float(fr.index()) / 1.1234;
                    *dr = double(fr.index()) / 1.1234;
                    *f1 = float(fr.index()) / 1.1234;
                    *d1 = double(fr.index()) / 1.1234;
                }
                convert<Tag_>(dvf2, dvd1);
                convert<Tag_>(dvd2, dvf1);
                dvf1.lock(lm_read_only, Tag_::memory_value);
                dvf1.unlock(lm_read_only);
                dvd1.lock(lm_read_only, Tag_::memory_value);
                dvd1.unlock(lm_read_only);
                dvf3.lock(lm_read_only, Tag_::memory_value);
                dvf3.unlock(lm_read_only);
                dvd3.lock(lm_read_only, Tag_::memory_value);
                dvd3.unlock(lm_read_only);
                convert<Tag_>(dvf3, dvd1);
                convert<Tag_>(dvd3, dvf1);
                convert<Tag_>(dvdf, dvf3);
                convert<Tag_>(dvff, dvd3);

                dvf2.lock(lm_read_only);
                dvd2.lock(lm_read_only);
                dvf2.unlock(lm_read_only);
                dvd2.unlock(lm_read_only);
                dvf3.lock(lm_read_only);
                dvd3.lock(lm_read_only);
                dvf3.unlock(lm_read_only);
                dvd3.unlock(lm_read_only);
                TEST_CHECK_EQUAL(dvf2, dvfr);
                TEST_CHECK_EQUAL(dvf3, dvfr);
                dvff.lock(lm_read_only);
                dvff.unlock(lm_read_only);
                dvdf.lock(lm_read_only);
                dvdf.unlock(lm_read_only);
                TEST_CHECK_EQUAL(dvff, dvfr);
                for (DenseVector<double>::ConstElementIterator i(dvd2.begin_elements()), r(dvdr.begin_elements()), i_end(dvd2.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
                for (DenseVector<double>::ConstElementIterator i(dvd3.begin_elements()), r(dvdr.begin_elements()), i_end(dvd3.end_elements())
                        ; i != i_end ; ++i, ++r)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
                }
                for (DenseVector<double>::ConstElementIterator i(dvdf.begin_elements()), r(dvdr.begin_elements()), i_end(dvdf.end_elements())
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
};
DenseVectorConvertTest<tags::CPU> dense_vector_convert_test;
#ifdef HONEI_SSE
DenseVectorConvertTest<tags::CPU::SSE> sse_dense_vector_convert_test;
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
DenseVectorConvertTest<tags::GPU::CUDA> cuda_dense_vector_convert_test;
#endif
#endif

template <typename Tag_>
class DenseVectorConvertQuickTest :
    public QuickTest
{
    public:
        DenseVectorConvertQuickTest() :
            QuickTest("dense_vector_convert_quick_test")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(4711);
            DenseVector<float> dvf1(size), dvf2(size, 4711), dvfr(size), dvf3(size);
            DenseVector<double> dvd1(size), dvd2(size, 4711), dvdr(size), dvd3(size);
            DenseVector<float>::ElementIterator fr(dvfr.begin_elements());
            DenseVector<double>::ElementIterator dr(dvdr.begin_elements());
            DenseVector<float>::ElementIterator f1(dvf1.begin_elements());
            DenseVector<double>::ElementIterator d1(dvd1.begin_elements());
            for (DenseVector<float>::ElementIterator i_end(dvf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
            {
                *fr = float(fr.index()) / 1.1234;
                *dr = double(fr.index()) / 1.1234;
                *f1 = float(fr.index()) / 1.1234;
                *d1 = double(fr.index()) / 1.1234;
            }
            convert<Tag_>(dvf2, dvd1);
            convert<Tag_>(dvd2, dvf1);

            dvf1.lock(lm_read_only, Tag_::memory_value);
            dvf1.unlock(lm_read_only);
            dvd1.lock(lm_read_only, Tag_::memory_value);
            dvd1.unlock(lm_read_only);
            dvf3.lock(lm_read_only, Tag_::memory_value);
            dvf3.unlock(lm_read_only);
            dvd3.lock(lm_read_only, Tag_::memory_value);
            dvd3.unlock(lm_read_only);
            convert<Tag_>(dvf3, dvd1);
            convert<Tag_>(dvd3, dvf1);

            dvf2.lock(lm_read_only);
            dvd2.lock(lm_read_only);
            dvf2.unlock(lm_read_only);
            dvd2.unlock(lm_read_only);
            dvf3.lock(lm_read_only);
            dvd3.lock(lm_read_only);
            dvf3.unlock(lm_read_only);
            dvd3.unlock(lm_read_only);
            TEST_CHECK_EQUAL(dvf2, dvfr);
            TEST_CHECK_EQUAL(dvf3, dvfr);
            for (DenseVector<double>::ConstElementIterator i(dvd2.begin_elements()), r(dvdr.begin_elements()), i_end(dvd2.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            for (DenseVector<double>::ConstElementIterator i(dvd3.begin_elements()), r(dvdr.begin_elements()), i_end(dvd3.end_elements())
                    ; i != i_end ; ++i, ++r)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *r, sqrt(std::numeric_limits<float>::epsilon()));
            }
            DenseVector<float> dvf01(3), dvf02(3);
            DenseVector<double> dvd01(2), dvd02(4);
            TEST_CHECK_THROWS(convert(dvf01, dvd01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(convert(dvd02, dvf02), VectorSizeDoesNotMatch);
        }
};
DenseVectorConvertQuickTest<tags::CPU> dense_vector_convert_quick_test;
#ifdef HONEI_SSE
DenseVectorConvertQuickTest<tags::CPU::SSE> sse_dense_vector_convert_quick_test;
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
DenseVectorConvertQuickTest<tags::GPU::CUDA> cuda_dense_vector_convert_quick_test;
#endif
#endif

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
                        DenseVector<float>::ElementIterator fr(dvfrr.begin_elements());
                        DenseVector<double>::ElementIterator dr(dvdrr.begin_elements());
                        DenseVector<float>::ElementIterator f1(dvf1r.begin_elements());
                        DenseVector<double>::ElementIterator d1(dvd1r.begin_elements());
                        for (DenseVector<float>::ElementIterator i_end(dvf1r.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                        {
                            *fr = float(fr.index()) / 1.1234;
                            *dr = double(fr.index()) / 1.1234;
                            *f1 = float(fr.index()) / 1.1234;
                            *d1 = double(fr.index()) / 1.1234;
                        }
                        convert(dvf2r, dvd1r);
                        convert(dvd2r, dvf1r);
                        TEST_CHECK_EQUAL(dvf2r, dvfrr);
                        for (DenseVector<double>::ConstElementIterator i(dvd2r.begin_elements()), r(dvdrr.begin_elements()), i_end(dvd2r.end_elements())
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
                    DenseVector<float>::ElementIterator fr(dvfrr.begin_elements());
                    DenseVector<double>::ElementIterator dr(dvdrr.begin_elements());
                    DenseVector<float>::ElementIterator f1(dvf1r.begin_elements());
                    DenseVector<double>::ElementIterator d1(dvd1r.begin_elements());
                    for (DenseVector<float>::ElementIterator i_end(dvf1r.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
                    {
                        *fr = float(fr.index()) / 1.1234;
                        *dr = double(fr.index()) / 1.1234;
                        *f1 = float(fr.index()) / 1.1234;
                        *d1 = double(fr.index()) / 1.1234;
                    }
                    convert(dvf2r, dvd1r);
                    convert(dvd2r, dvf1r);
                    TEST_CHECK_EQUAL(dvf2r, dvfrr);
                    for (DenseVector<double>::ConstElementIterator i(dvd2r.begin_elements()), r(dvdrr.begin_elements()), i_end(dvd2r.end_elements())
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
                SparseMatrix<float>::ElementIterator fr(smfr.begin_elements());
                SparseMatrix<double>::ElementIterator dr(smdr.begin_elements());
                SparseMatrix<float>::ElementIterator f1(smf1.begin_elements());
                SparseMatrix<double>::ElementIterator d1(smd1.begin_elements());
                for (SparseMatrix<float>::ElementIterator i_end(smf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
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
                for (SparseMatrix<double>::ConstElementIterator i(smd2.begin_elements()), r(smdr.begin_elements()), i_end(smd2.end_elements())
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
            SparseMatrix<float>::ElementIterator fr(smfr.begin_elements());
            SparseMatrix<double>::ElementIterator dr(smdr.begin_elements());
            SparseMatrix<float>::ElementIterator f1(smf1.begin_elements());
            SparseMatrix<double>::ElementIterator d1(smd1.begin_elements());
            for (SparseMatrix<float>::ElementIterator i_end(smf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
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
            for (SparseMatrix<double>::ConstElementIterator i(smd2.begin_elements()), r(smdr.begin_elements()), i_end(smd2.end_elements())
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
                SparseVector<float>::ElementIterator fr(svfr.begin_elements());
                SparseVector<double>::ElementIterator dr(svdr.begin_elements());
                SparseVector<float>::ElementIterator f1(svf1.begin_elements());
                SparseVector<double>::ElementIterator d1(svd1.begin_elements());
                for (SparseVector<float>::ElementIterator i_end(svf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
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
                for (SparseVector<double>::ConstElementIterator i(svd2.begin_elements()), r(svdr.begin_elements()), i_end(svd2.end_elements())
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
            SparseVector<float>::ElementIterator fr(svfr.begin_elements());
            SparseVector<double>::ElementIterator dr(svdr.begin_elements());
            SparseVector<float>::ElementIterator f1(svf1.begin_elements());
            SparseVector<double>::ElementIterator d1(svd1.begin_elements());
            for (SparseVector<float>::ElementIterator i_end(svf1.end_elements()) ; f1 != i_end ; ++fr, ++dr, ++f1, ++d1)
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
            for (SparseVector<double>::ConstElementIterator i(svd2.begin_elements()), r(svdr.begin_elements()), i_end(svd2.end_elements())
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

template <typename Tag_, typename DT_>
class DenseVectorContinuousBaseFillQuickTest :
    public QuickTest
{
    public:
        DenseVectorContinuousBaseFillQuickTest(const std::string & type) :
            QuickTest("dense_vector_continuous_base_fill_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(4711);
            DenseVector<DT_> dv(2 * size);
            DenseVectorRange<DT_> dvr1(dv, size / 3, 0);
            DenseVectorRange<DT_> dvr2(dv, size, size / 3);
            DenseVectorRange<DT_> dvr3(dv, size - 2 * size / 3, size + size / 3);

            fill<Tag_>(dv, DT_(8.472));

            for (typename DenseVector<DT_>::ConstElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(8.472), std::numeric_limits<DT_>::epsilon());
            }

            fill<Tag_>(dvr2, DT_(75.633));

            for (typename DenseVectorRange<DT_>::ConstElementIterator i(dvr1.begin_elements()), i_end(dvr1.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(8.472), std::numeric_limits<DT_>::epsilon());
            }

            for (typename DenseVectorRange<DT_>::ConstElementIterator i(dvr2.begin_elements()), i_end(dvr2.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(75.633), std::numeric_limits<DT_>::epsilon());
            }

            for (typename DenseVectorRange<DT_>::ConstElementIterator i(dvr3.begin_elements()), i_end(dvr3.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(8.472), std::numeric_limits<DT_>::epsilon());
            }
        }
};

DenseVectorContinuousBaseFillQuickTest<tags::CPU, float> dense_vector_continuous_base_fill_quick_test_float("float");
DenseVectorContinuousBaseFillQuickTest<tags::CPU, double> dense_vector_continuous_base_fill_quick_test_double("double");

template <typename Tag_, typename DT_>
class DenseMatrixFillQuickTest :
    public QuickTest
{
    public:
        DenseMatrixFillQuickTest(const std::string & type) :
            QuickTest("dense_matrix_fill_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long rows(47), columns(11);
            DenseMatrix<DT_> dm(rows, columns);

            fill<Tag_>(dm, DT_(8.472));

            for (typename DenseMatrix<DT_>::ConstElementIterator i(dm.begin_elements()), i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(8.472), std::numeric_limits<DT_>::epsilon());
            }
        }
};

DenseMatrixFillQuickTest<tags::CPU, float> dense_matrix_fill_quick_test_float("float");
DenseMatrixFillQuickTest<tags::CPU, double> dense_matrix_fill_quick_test_double("double");

template <typename Tag_, typename DT_>
class CudaDenseMatrixFillQuickTest :
    public QuickTest
{
    public:
        CudaDenseMatrixFillQuickTest(const std::string & type) :
            QuickTest("cuda_dense_matrix_fill_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long rows(47), columns(11);
            DenseMatrix<DT_> dm(rows, columns, DT_(25.3));

            fill<Tag_>(dm, DT_(0));

            dm.lock(lm_read_only);
            for (typename DenseMatrix<DT_>::ConstElementIterator i(dm.begin_elements()), i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DT_(0), std::numeric_limits<DT_>::epsilon());
            }
            dm.unlock(lm_read_only);
        }
};
#ifdef HONEI_CUDA
CudaDenseMatrixFillQuickTest<tags::GPU::CUDA, float> cuda_dense_matrix_fill_quick_test_float("float");
#endif

template <typename DT_>
class IteratorFillQuickTest :
    public QuickTest
{
    public:
        IteratorFillQuickTest(const std::string & type) :
            QuickTest("iterator_fill_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(123);
            DenseVector<DT_> dv1(size, 1.2), dv2(size, DT_(1.2));
            DenseVector<DT_> dvs(dv1, size / 5, size / 3, 2);

            fill(dvs.begin_elements(), dvs.end_elements(), DT_(3.4));

            for (typename DenseVector<DT_>::ElementIterator i(dv2.element_at(size / 3)), i_end(dv2.element_at(size / 3 + 2 * (size / 5))) ;
                    i != i_end ; i += 2)
            {
                *i = DT_(3.4);
            }

            for (typename DenseVector<DT_>::ConstElementIterator i(dv1.begin_elements()), i_end(dv1.begin_elements()), j(dv1.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DT_>::epsilon());
            }
        }
};
/// \todo Instanziieren

template <typename Tag_, typename DT_>
class DenseVectorCopyQuickTest :
    public QuickTest
{
    public:
        DenseVectorCopyQuickTest(const std::string & type) :
            QuickTest("dense_vector_copy_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(47);
            DenseVector<DT_> dv(size);
            for (unsigned long i(0) ; i < dv.size() ; ++i)
            {
                dv[i] = DT_(i % 10);
            }
            DenseVector<DT_> dv2(size, DT_(471));
            copy<Tag_>(dv, dv2);
            dv.lock(lm_read_only);
            dv2.lock(lm_read_only);
            TEST_CHECK_EQUAL(dv2, dv);
            dv.unlock(lm_read_only);
            dv2.unlock(lm_read_only);
            copy<Tag_>(dv2, dv);
            copy<Tag_>(dv, dv2);
            dv.lock(lm_read_only);
            dv2.lock(lm_read_only);
            TEST_CHECK_EQUAL(dv2, dv);
            dv.unlock(lm_read_only);
            dv2.unlock(lm_read_only);


            dv.lock(lm_write_only);
            dv2.lock(lm_write_only);
            for (unsigned long i(0) ; i < dv.size() ; ++i)
            {
                dv[i] = DT_(i % 10);
            }
            for (unsigned long i(0) ; i < dv.size() ; ++i)
            {
                dv2[i] = DT_(i);
            }
            dv.unlock(lm_write_only);
            dv2.unlock(lm_write_only);
            dv.lock(lm_read_only, Tag_::memory_value);
            dv2.lock(lm_read_only, Tag_::memory_value);
            dv.unlock(lm_read_only);
            dv2.unlock(lm_read_only);
            copy<Tag_>(dv, dv2);
            dv.lock(lm_read_only);
            dv2.lock(lm_read_only);
            TEST_CHECK_EQUAL(dv2, dv);
            dv.unlock(lm_read_only);
            dv2.unlock(lm_read_only);
        }
};
DenseVectorCopyQuickTest<tags::CPU, float> dense_vector_copy_quick_test_float("float");
DenseVectorCopyQuickTest<tags::CPU, double> dense_vector_copy_quick_test_double("double");
#ifdef HONEI_CUDA
DenseVectorCopyQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_copy_quick_test_float("float");
DenseVectorCopyQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_copy_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorCopyQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_copy_quick_test_double("double");
DenseVectorCopyQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_copy_quick_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorCopyQuickTest<tags::OpenCL::CPU, float> ocl_cpu_dense_vector_copy_quick_test_float("float");
DenseVectorCopyQuickTest<tags::OpenCL::CPU, double> ocl_cpu_dense_vector_copy_quick_test_double("double");
#ifdef HONEI_CUDA
DenseVectorCopyQuickTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_copy_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorCopyQuickTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_copy_quick_test_double("double");
#endif
#endif
#endif
