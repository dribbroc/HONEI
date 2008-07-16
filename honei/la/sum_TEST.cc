/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
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

#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/sum.hh>
#include <unittest/unittest.hh>
#include <honei/util/memory_arbiter.hh>

#include <limits>
#include <tr1/memory>
#include <cstdlib>

using namespace honei;
using namespace tests;

// Test cases for matrix operations

template <typename Tag_, typename DataType_>
class BandedMatrixDenseMatrixSumTest :
    public BaseTest
{
    private: 
        DataType_ getRandom() const
        {
            return DataType_(std::rand()) / DataType_(std::rand()-RAND_MAX / 2);
        }

    public:
        BandedMatrixDenseMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size);
                //randomize content
                for (typename DenseVector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()); i != i_end; ++i)
                {
                    *i =getRandom();
                }

                // multiple bands (with same vector, perhaps should be different)
                BandedMatrix<DataType_> bm1(size, dv1);
                bm1.insert_band( 2, dv1);
                bm1.insert_band(-3, dv1);
                bm1.insert_band( 4, dv1);
                bm1.insert_band(-5, dv1);
                bm1.insert_band(-(size - 1), dv1);  // Boundary element, lower left.

                // randomized content for dense matrices
                DenseMatrix<DataType_> dm2(size, size), dm3(size, size);
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm2.begin_elements()), i_end(dm2.end_elements()), j(dm3.begin_elements());
                    i != i_end; ++i, ++j)
                {
                    *i = *j = getRandom();
                } 

                // calculate reference result
                typename MutableMatrix<DataType_>::ElementIterator k(dm3.begin_elements());
                for (typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                    i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
                { 
                    if (*i != DataType_(0))
                     {
                        *k = DataType_(*i) + dm2[i.row()][i.column()];
                    }
                }

                Sum<Tag_>::value(dm2, bm1);
                TEST_CHECK_EQUAL(dm2, dm3);
            }

            BandedMatrix<DataType_> bm01(5);
            DenseMatrix<DataType_> dm02(6, 6);

            TEST_CHECK_THROWS(Sum<Tag_>::value(dm02, bm01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseMatrixSumTest<tags::CPU, float> banded_matrix_dense_matrix_sum_test_float("float");
BandedMatrixDenseMatrixSumTest<tags::CPU, double> banded_matrix_dense_matrix_sum_test_double("double");
BandedMatrixDenseMatrixSumTest<tags::CPU::MultiCore, float> mc_banded_matrix_dense_matrix_sum_test_float("MC float");
BandedMatrixDenseMatrixSumTest<tags::CPU::MultiCore, double> mc_banded_matrix_dense_matrix_sum_test_double("MC double");

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
            DenseVector<DataType_> dv1(size, DataType_(2));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseMatrix<DataType_> dm2(size, size, DataType_(1)), dm3(size, size, DataType_(1));

            typename MutableMatrix<DataType_>::ElementIterator k(dm3.begin_elements());
            for (typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k)
            {
                if (*i != DataType_(0))
                {
                    *k = DataType_(3);
                }
            }
            Sum<>::value(dm2, bm1);

            TEST_CHECK_EQUAL(dm2, dm3);

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
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
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
                typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements());
                for (typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                        i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k, ++j)
                {
                    if (*i != DataType_(0) && *j != DataType_(0))
                    {
                        *k = DataType_(3);
                    }
                    else if (*i != DataType_(0))
                    {
                        *k = DataType_(2);
                    }
                }
                Sum<>::value(sm2, bm1);

                TEST_CHECK_EQUAL(sm2, sm3);

                BandedMatrix<DataType_> bm01(5);
                SparseMatrix<DataType_> sm02(6, 6, 1), sm03(5, 6, 1);

                TEST_CHECK_THROWS(Sum<>::value(sm02, bm01), MatrixRowsDoNotMatch);
                TEST_CHECK_THROWS(Sum<>::value(sm03, bm01), MatrixIsNotSquare);
            }
        }
};
BandedMatrixSparseMatrixSumTest<float> banded_matrix_sparse_matrix_sum_test_float("float");
BandedMatrixSparseMatrixSumTest<double> banded_matrix_sparse_matrix_sum_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_matrix_sum_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(20); 
            DenseVector<DataType_> dv1(size, DataType_(2));
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
            typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements());
            for (typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements()),
                    i_end(bm1.end_elements()) ; i != i_end ; ++i, ++k, ++j)
            {
                if (*i != DataType_(0) && *j != DataType_(0))
                {
                    *k = DataType_(3);
                }
                else if (*i != DataType_(0))
                {
                    *k = DataType_(2);
                }
            }
            Sum<>::value(sm2, bm1);

            TEST_CHECK_EQUAL(sm2, sm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(Sum<>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<>::value(sm03, bm01), MatrixIsNotSquare);
        }
};
BandedMatrixSparseMatrixSumQuickTest<float> banded_matrix_sparse_matrix_sum_quick_test_float("float");
BandedMatrixSparseMatrixSumQuickTest<double> banded_matrix_sparse_matrix_sum_quick_test_double("double");

template <typename Tag_, typename DataType_>
class BandedMatrixSumTest :
    public BaseTest
{
    public:
        BandedMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
                Sum<Tag_>::value(bm1, bm2);
                for (typename BandedMatrix<DataType_>::ConstBandIterator ce(bm1.begin_bands()),
                        ce_end(bm1.end_bands()) ; ce != ce_end ; ++ce)
                {
                    if (ce.index() != size - 1 )
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                                i_end(ce->end_elements()) ; i < i_end ; ++i)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                        }
                        }
                        else
                        {
                            for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                                    i_end(ce->end_elements()) ; i < i_end ; ++i)
                            {
                                TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                            }
                        }
                    }
                }

                BandedMatrix<DataType_> bm01(5), bm02(6);

                TEST_CHECK_THROWS(Sum<Tag_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixSumTest<tags::CPU, float> banded_matrix_sum_test_float("float");
BandedMatrixSumTest<tags::CPU, double> banded_matrix_sum_test_double("double");
BandedMatrixSumTest<tags::CPU::MultiCore, float> mc_banded_matrix_sum_test_float("MC float");
BandedMatrixSumTest<tags::CPU::MultiCore, double> mc_banded_matrix_sum_test_double("MC double");


template <typename Tag_, typename DataType_>
class BandedMatrixSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv2(size, DataType_(3));
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
            Sum<Tag_>::value(bm1, bm2);
            for (typename BandedMatrix<DataType_>::ConstBandIterator ce(bm1.begin_bands()), ce_end(bm1.end_bands()) ;
                    ce != ce_end ; ++ce)
            {
                if (ce.index() != size - 1 )
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                            i_end(ce->end_elements()) ; i < i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    }
                }
                else
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                        i_end(ce->end_elements()) ; i < i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                    }
                }
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(Sum<Tag_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixSumQuickTest<tags::CPU, float> banded_matrix_sum_quick_test_float("float");
BandedMatrixSumQuickTest<tags::CPU, double> banded_matrix_sum_quick_test_double("double");
BandedMatrixSumQuickTest<tags::CPU::MultiCore, float> mc_banded_matrix_sum_quick_test_float("MC float");
BandedMatrixSumQuickTest<tags::CPU::MultiCore, double> mc_banded_matrix_sum_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseMatrixSparseMatrixSumTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixSumTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(0)),
                        dm3(size+1, size, DataType_(0));
                SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1);
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
                Sum<Tag_>::value(dm1, sm2);

                TEST_CHECK_EQUAL(dm1, dm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DataType_> dm03(6, 5);

            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixSumTest<tags::CPU, float> dense_matrix_sparse_matrix_sum_test_float("float");
DenseMatrixSparseMatrixSumTest<tags::CPU, double> dense_matrix_sparse_matrix_sum_test_double("double");
DenseMatrixSparseMatrixSumTest<tags::CPU::MultiCore, float> mc_dense_matrix_sparse_matrix_sum_test_float("MC float");
DenseMatrixSparseMatrixSumTest<tags::CPU::MultiCore, double> mc_dense_matrix_sparse_matrix_sum_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseMatrixSparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(0)),
                    dm3(size+1, size, DataType_(0));
            SparseMatrix<DataType_> sm2(size+1, size, size / 7 + 1);
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
            Sum<Tag_>::value(dm1, sm2);

            TEST_CHECK_EQUAL(dm1, dm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1);
            DenseMatrix<DataType_> dm03(6, 5);

            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, sm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSparseMatrixSumQuickTest<tags::CPU, float> dense_matrix_sparse_matrix_sum_quick_test_float("float");
DenseMatrixSparseMatrixSumQuickTest<tags::CPU, double> dense_matrix_sparse_matrix_sum_quick_test_double("double");
DenseMatrixSparseMatrixSumQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_sparse_matrix_sum_quick_test_float("MC float");
DenseMatrixSparseMatrixSumQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_sparse_matrix_sum_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseMatrixSumTest :
    public BaseTest
{
    public:
        DenseMatrixSumTest(const std::string & type) :
            BaseTest("dense_matrix_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(-3)),
                    dm3(size+1, size, DataType_(-1));

                Sum<Tag_>::value(dm1, dm2);

                TEST(dm1.lock(lm_read_only), TEST_CHECK_EQUAL(dm1, dm3), dm1.unlock(lm_read_only));

                DenseMatrix<DataType_> dm4(size + 1, size, DataType_(117));
                dm2.lock(lm_read_and_write);

                for (typename MutableMatrix<DataType_>::ElementIterator i(dm4.begin_elements()), i_end(dm4.end_elements()),
                        si(dm2.begin_elements()), ri(dm3.begin_elements()) ; i != i_end ; ++i, ++si, ++ri)
                {
                    *i  = *i + i.row() + i.column();
                    *si = *si + si.row() - si.column();
                    *ri = DataType_(117) + DataType_(2) * ri.row() - DataType_(3);
                }

                dm2.unlock(lm_read_and_write);
                Sum<Tag_>::value(dm4, dm2);
                TEST(dm4.lock(lm_read_only), TEST_CHECK_EQUAL(dm4, dm3), dm4.unlock(lm_read_only));
            }

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(6, 5);

            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSumTest<tags::CPU, float> dense_matrix_sum_test_float("float");
DenseMatrixSumTest<tags::CPU, double> dense_matrix_sum_test_double("double");
DenseMatrixSumTest<tags::CPU::MultiCore, float> mc_dense_matrix_sum_test_float("MC float");
DenseMatrixSumTest<tags::CPU::MultiCore, double> mc_dense_matrix_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixSumTest<tags::CPU::SSE, float> sse_dense_matrix_sum_test_float("SSE float");
DenseMatrixSumTest<tags::CPU::SSE, double> sse_dense_matrix_sum_test_double("SSE double");
DenseMatrixSumTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_sum_test_float("MC SSE float");
DenseMatrixSumTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseMatrixSumTest<tags::GPU::CUDA, float> cuda_dense_matrix_sum_test_float("float");
#endif
#ifdef HONEI_CELL
DenseMatrixSumTest<tags::Cell, float> cell_dense_matrix_sum_test_float("Cell float");
DenseMatrixSumTest<tags::Cell, double> cell_dense_matrix_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(-3)),
                dm3(size+1, size, DataType_(-1));
            Sum<Tag_>::value(dm1, dm2);

            TEST(dm1.lock(lm_read_only), TEST_CHECK_EQUAL(dm1, dm3), dm1.unlock(lm_read_only));

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(6, 5);

            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixSumQuickTest<tags::CPU, float> dense_matrix_sum_quick_test_float("float");
DenseMatrixSumQuickTest<tags::CPU, double> dense_matrix_sum_quick_test_double("double");
DenseMatrixSumQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_sum_quick_test_float("MC float");
DenseMatrixSumQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixSumQuickTest<tags::CPU::SSE, float> sse_dense_matrix_sum_quick_test_float("SSE float");
DenseMatrixSumQuickTest<tags::CPU::SSE, double> sse_dense_matrix_sum_quick_test_double("SSE double");
DenseMatrixSumQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_sum_quick_test_float("MC SSE float");
DenseMatrixSumQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseMatrixSumQuickTest<tags::GPU::CUDA, float> cuda_dense_matrix_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseMatrixSumQuickTest<tags::Cell, float> cell_dense_matrix_sum_quick_test_float("Cell float");
DenseMatrixSumQuickTest<tags::Cell, double> cell_dense_matrix_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixSumTest :
    public BaseTest
{
    public:
        SparseMatrixSumTest(const std::string & type) :
            BaseTest("sparse_matrix_sum_test<" + type + ">")
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
                Sum<Tag_>::value(sm1, sm2);

                TEST_CHECK_EQUAL(sm1, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(Sum<Tag_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixSumTest<tags::CPU, float> sparse_matrix_sum_test_float("float");
SparseMatrixSumTest<tags::CPU, double> sparse_matrix_sum_test_double("double");
SparseMatrixSumTest<tags::CPU::MultiCore, float> mc_sparse_matrix_sum_test_float("MC float");
SparseMatrixSumTest<tags::CPU::MultiCore, double> mc_sparse_matrix_sum_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        SparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_sum_quick_test<" + type + ">")
    {
            register_tag(Tag_::name);
    }

        virtual void run() const
        {
            unsigned long size (11);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                sm2(size+1, size, size / 7 + 1), sm3(size+1, size, size / 8 + 1 );
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
            Sum<Tag_>::value(sm1, sm2);

            TEST_CHECK_EQUAL(sm1, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(Sum<Tag_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(Sum<Tag_>::value(sm03, sm02), MatrixColumnsDoNotMatch);
        }
};
SparseMatrixSumQuickTest<tags::CPU, float> sparse_matrix_sum_quick_test_float("float");
SparseMatrixSumQuickTest<tags::CPU, double> sparse_matrix_sum_quick_test_double("double");
SparseMatrixSumQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_sum_quick_test_float("MC float");
SparseMatrixSumQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_sum_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class ScalarDenseMatrixSumTest :
    public BaseTest
{
    public:
        ScalarDenseMatrixSumTest(const std::string & type) :
            BaseTest("scalar_dense_matrix_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(5));
                Sum<Tag_>::value(dm1, DataType_(3));

                TEST_CHECK_EQUAL(dm1, dm2);

                DenseMatrix<DataType_> dm3(size + 1, size, DataType_(0));
                DenseMatrix<DataType_> dm4(size + 1, size, DataType_(0)),
                                       dm5(size + 1, size, DataType_(0)),
                                       dm6(size + 1, size, DataType_(0));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm3.begin_elements()),
                        i_end(dm3.end_elements()), li(dm4.begin_elements()),
                        mi(dm5.begin_elements()), ni(dm6.begin_elements())
                        ; i != i_end ; ++i, ++li, ++mi, ++ni)
                {
                    *i = i.row() + i.column() - DataType_(size / 2);
                    *li = *i;
                    *mi = *i - 7;
                    *ni = *i + 16501;
                }
                Sum<Tag_>::value(dm3, DataType_(0));
                TEST_CHECK_EQUAL(dm3, dm4);
                Sum<Tag_>::value(dm3, DataType_(-7));
                TEST_CHECK_EQUAL(dm3, dm5);
                Sum<Tag_>::value(dm3, DataType_(16508));
                TEST_CHECK_EQUAL(dm3, dm6);
            }
        }
};
ScalarDenseMatrixSumTest<tags::CPU, float> scalar_dense_matrix_sum_test_float("float");
ScalarDenseMatrixSumTest<tags::CPU, double> scalar_dense_matrix_sum_test_double("double");
#ifdef HONEI_SSE
ScalarDenseMatrixSumTest<tags::CPU::SSE, float> sse_scalar_dense_matrix_sum_test_float("SSE float");
ScalarDenseMatrixSumTest<tags::CPU::SSE, double> sse_scalar_dense_matrix_sum_test_double("SSE double");
ScalarDenseMatrixSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_matrix_sum_test_float("MC SSE float");
ScalarDenseMatrixSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_matrix_sum_test_double("MC SSE double");
#endif
ScalarDenseMatrixSumTest<tags::CPU::MultiCore, float> mc_scalar_dense_matrix_sum_test_float("MC float");
ScalarDenseMatrixSumTest<tags::CPU::MultiCore, double> mc_scalar_dense_matrix_sum_test_double("MC double");
#ifdef HONEI_CELL
ScalarDenseMatrixSumTest<tags::Cell, float> cell_scalar_dense_matrix_sum_test_float("Cell float");
ScalarDenseMatrixSumTest<tags::Cell, double> cell_scalar_dense_matrix_sum_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class ScalarDenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_matrix_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm(size+1, size, DataType_(2));
            Sum<Tag_>::value(dm, DataType_(3));

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
ScalarDenseMatrixSumQuickTest<tags::CPU, float> scalar_dense_matrix_sum_quick_test_float("float");
ScalarDenseMatrixSumQuickTest<tags::CPU, double> scalar_dense_matrix_sum_quick_test_double("double");
#ifdef HONEI_SSE
ScalarDenseMatrixSumQuickTest<tags::CPU::SSE, float> sse_scalar_dense_matrix_sum_quick_test_float("SSE float");
ScalarDenseMatrixSumQuickTest<tags::CPU::SSE, double> sse_scalar_dense_matrix_sum_quick_test_double("SSE double");
ScalarDenseMatrixSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_matrix_sum_quick_test_float("MC SSE float");
ScalarDenseMatrixSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_matrix_sum_quick_test_double("MCSSE double");
#endif
ScalarDenseMatrixSumQuickTest<tags::CPU::MultiCore, float> mc_scalar_dense_matrix_sum_quick_test_float("MC float");
ScalarDenseMatrixSumQuickTest<tags::CPU::MultiCore, double> mc_scalar_dense_matrix_sum_quick_test_double("MC double");
#ifdef HONEI_CELL
ScalarDenseMatrixSumQuickTest<tags::Cell, float> cell_scalar_dense_matrix_sum_quick_test_float("Cell float");
ScalarDenseMatrixSumQuickTest<tags::Cell, double> cell_scalar_dense_matrix_sum_quick_test_double("Cell double");
#endif

// Test cases for vector operations

template <typename Tag_, typename DataType_>
class ScalarDenseVectorSumTest :
    public BaseTest
{
    public:
        ScalarDenseVectorSumTest(const std::string & type) :
            BaseTest("scalar_dense_vector_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(9));
                Sum<Tag_>::value(dv1, DataType_(7));

                TEST_CHECK_EQUAL(dv1, dv2);
            }
        }
};
ScalarDenseVectorSumTest<tags::CPU, float>  scalar_dense_vector_sum_test_float("float");
ScalarDenseVectorSumTest<tags::CPU, double> scalar_dense_vector_sum_test_double("double");
#ifdef HONEI_SSE
ScalarDenseVectorSumTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_sum_test_float("SSE float");
ScalarDenseVectorSumTest<tags::CPU::SSE, double> sse_scalar_dense_vector_sum_test_double("SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorSumTest<tags::Cell, float> cell_scalar_dense_vector_sum_test_float("Cell float");
ScalarDenseVectorSumTest<tags::Cell, double> cell_scalar_dense_vector_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class ScalarDenseVectorSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_vector_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size = 18;
            DenseVector<DataType_> dv1(size, DataType_(4)), dv2(size, DataType_(9));
            Sum<Tag_>::value(dv1, DataType_(5));

            TEST_CHECK_EQUAL(dv1, dv2);
        }
};
ScalarDenseVectorSumQuickTest<tags::CPU, float>  scalar_dense_vector_sum_quick_test_float("float");
ScalarDenseVectorSumQuickTest<tags::CPU, double> scalar_dense_vector_sum_quick_test_double("double");
#ifdef HONEI_SSE
ScalarDenseVectorSumQuickTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_sum_quick_test_float("SSE float");
ScalarDenseVectorSumQuickTest<tags::CPU::SSE, double> sse_scalar_dense_vector_sum_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorSumQuickTest<tags::Cell, float>  cell_scalar_dense_vector_sum_quick_test_float("Cell float");
ScalarDenseVectorSumQuickTest<tags::Cell, double>  cell_scalar_dense_vector_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class ScalarDenseVectorRangeSumTest :
    public BaseTest
{
    public:
        ScalarDenseVectorRangeSumTest(const std::string & type) :
            BaseTest("scalar_dense_vector_range_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                for (unsigned i(0) ; i < 4 ; i++)
                {
                    DenseVector<DataType_> dv1(size+4, DataType_(2)), dv2(size+4, DataType_(9));
                    DenseVectorRange<DataType_> dvr1(dv1, size, i), dvr2(dv2, size, i);

                    Sum<Tag_>::value(dvr1, DataType_(7));
                    TEST_CHECK_EQUAL(dvr1, dvr2);
                }
            }
        }
};
ScalarDenseVectorRangeSumTest<tags::CPU, float>  scalar_dense_vector_range_sum_test_float("float");
ScalarDenseVectorRangeSumTest<tags::CPU, double> scalar_dense_vector_range_sum_test_double("double");
#ifdef HONEI_SSE
ScalarDenseVectorRangeSumTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_range_sum_test_float("SSE float");
ScalarDenseVectorRangeSumTest<tags::CPU::SSE, double> sse_scalar_dense_vector_range_sum_test_double("SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorRangeSumTest<tags::Cell, float> cell_scalar_dense_vector_range_sum_test_float("Cell float");
ScalarDenseVectorRangeSumTest<tags::Cell, double> cell_scalar_dense_vector_range_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class ScalarDenseVectorRangeSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorRangeSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_vector_range_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size = 18;
            for (unsigned i(0) ; i < 4 ; i++)
            {
                DenseVector<DataType_> dv1(size+4, DataType_(2)), dv2(size+4, DataType_(9));
                DenseVectorRange<DataType_> dvr1(dv1, size, i), dvr2(dv2, size, i);

                Sum<Tag_>::value(dvr1, DataType_(7));
                TEST_CHECK_EQUAL(dvr1, dvr2);
            }
        }
};
ScalarDenseVectorRangeSumQuickTest<tags::CPU, float>  scalar_dense_vector_range_sum_quick_test_float("float");
ScalarDenseVectorRangeSumQuickTest<tags::CPU, double> scalar_dense_vector_range_sum_quick_test_double("double");
#ifdef HONEI_SSE
ScalarDenseVectorRangeSumQuickTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_range_sum_quick_test_float("SSE float");
ScalarDenseVectorRangeSumQuickTest<tags::CPU::SSE, double> sse_scalar_dense_vector_range_sum_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorRangeSumQuickTest<tags::Cell, float>  cell_scalar_dense_vector_range_sum_quick_test_float("Cell float");
ScalarDenseVectorRangeSumQuickTest<tags::Cell, double>  cell_scalar_dense_vector_range_sum_quick_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class DenseVectorRangeSumTest :
    public BaseTest
{
    public:
        DenseVectorRangeSumTest(const std::string & type) :
            BaseTest("dense_vector_range_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(40) ; size < (1 << 15) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size), dv2(size);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    *i = (i.index() + 1) / 1.23456789;
                    *j = 2 * *i;
                }

                DenseVectorRange<DataType_> dvr1 (dv1, size / 2, 5);
                DenseVectorRange<DataType_> dvr2 (dv2, size / 2, 5);
                Sum<Tag_>::value(dvr1, dvr2);

                TEST(dvr1.lock(lm_read_only),
                        for (typename Vector<DataType_>::ConstElementIterator i(dvr1.begin_elements()),
                            i_end(dvr1.end_elements()) ; i != i_end ; ++i)
                        {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_((i.index() + 5 + 1) * 3 / 1.23456789),
                            std::numeric_limits<DataType_>::epsilon() * 2 * (i.index() + 5 + 1) * 3 / 1.23456789);
                        },
                        dvr1.unlock(lm_read_only));
                DenseVectorRange<DataType_> dvr3 (dv1, size / 6, (size / 2) + 2);
                DenseVectorRange<DataType_> dvr4 (dv2, size / 6, size / 2);
                Sum<Tag_>::value(dvr3, dvr4);
            }

            DenseVector<DataType_> dv01(6, DataType_(1)), dv02(4, DataType_(1));
            DenseVectorRange<DataType_> dvr01(dv01, 4, 2), dvr02(dv02, 3 , 1);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dvr01, dvr02), VectorSizeDoesNotMatch);
        }
};
DenseVectorRangeSumTest<tags::CPU, float> dense_vector_range_sum_test_float("float");
DenseVectorRangeSumTest<tags::CPU, double> dense_vector_range_sum_test_double("double");
DenseVectorRangeSumTest<tags::CPU::MultiCore, float> mc_dense_vector_range_sum_test_float("MC float");
DenseVectorRangeSumTest<tags::CPU::MultiCore, double> mc_dense_vector_range_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeSumTest<tags::CPU::SSE, float> sse_dense_vector_range_sum_test_float("SSE float");
DenseVectorRangeSumTest<tags::CPU::SSE, double> sse_dense_vector_range_sum_test_double("SSE double");
DenseVectorRangeSumTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_sum_test_float("MC SSE float");
DenseVectorRangeSumTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeSumTest<tags::GPU::CUDA, float> cuda_dense_vector_range_sum_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeSumTest<tags::Cell, float> cell_dense_vector_range_sum_test_float("Cell float");
DenseVectorRangeSumTest<tags::Cell, double> cell_dense_vector_range_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(40);
            DenseVector<DataType_> dv1(size), dv2(size);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                *i = (i.index() + 1) / 1.23456789;
                *j = 2 * *i;
            }

            DenseVectorRange<DataType_> dvr1 (dv1, 11, 5);
            DenseVectorRange<DataType_> dvr2 (dv2, 11, 5);
            Sum<Tag_>::value(dvr1, dvr2);

            TEST(dvr1.lock(lm_read_only),
                    for (typename Vector<DataType_>::ConstElementIterator i(dvr1.begin_elements()),
                        i_end(dvr1.end_elements()) ; i != i_end ; ++i)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_((i.index() + 5 + 1) * 3 / 1.23456789),
                        std::numeric_limits<DataType_>::epsilon() * 2 * (i.index() + 5 + 1) * 3 / 1.23456789);
                    },
                    dvr1.unlock(lm_read_only));

            DenseVector<DataType_> dv01(6, DataType_(1)), dv02(4, DataType_(1));
            DenseVectorRange<DataType_> dvr01(dv01, 4, 2), dvr02(dv02, 3 , 1);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dvr01, dvr02), VectorSizeDoesNotMatch);
        }
};
DenseVectorRangeSumQuickTest<tags::CPU, float> dense_vector_range_sum_quick_test_float("float");
DenseVectorRangeSumQuickTest<tags::CPU, double> dense_vector_range_sum_quick_test_double("double");
DenseVectorRangeSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_range_sum_quick_test_float("MC float");
DenseVectorRangeSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_sum_quick_test_float("SSE float");
DenseVectorRangeSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_sum_quick_test_double("SSE double");
DenseVectorRangeSumQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_sum_quick_test_float("MC SSE float");
DenseVectorRangeSumQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_range_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeSumQuickTest<tags::Cell, float> cell_dense_vector_range_sum_quick_test_float("Cell float");
DenseVectorRangeSumQuickTest<tags::Cell, double> cell_dense_vector_range_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 15) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size + 3), dv2(size + 3);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    *i = (i.index() + 1) / 1.23456789;
                    *j = 2 * *i;
                }

                Sum<Tag_>::value(dv1, dv2);

                dv1.lock(lm_read_only);
                for (typename Vector<DataType_>::ConstElementIterator i(dv1.begin_elements()),
                        i_end(dv1.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_((i.index() + 1) * 3 / 1.23456789),
                            std::numeric_limits<DataType_>::epsilon() * 2 * (i.index() + 1) * 3 / 1.23456789);
                }
                dv1.unlock(lm_read_only);

                Sum<Tag_>::value(dv1, dv2);

                TEST(dv1.lock(lm_read_only),
                        for (typename Vector<DataType_>::ConstElementIterator i(dv1.begin_elements()),
                            i_end(dv1.end_elements()) ; i != i_end ; ++i)
                        {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_((i.index() + 1) * 5 / 1.23456789),
                            std::numeric_limits<DataType_>::epsilon() * 2 * (i.index() + 1) * 5 / 1.23456789);
                        },
                        dv1.unlock(lm_read_only));

                DenseVector<DataType_> dv01(6, DataType_(1)), dv02(4, DataType_(1));
                TEST_CHECK_THROWS(Sum<Tag_>::value(dv01, dv02), VectorSizeDoesNotMatch);
            }
        }
};
DenseVectorSumTest<tags::CPU, float> dense_vector_sum_test_float("float");
DenseVectorSumTest<tags::CPU, double> dense_vector_sum_test_double("double");
DenseVectorSumTest<tags::CPU::MultiCore, float> mc_dense_vector_sum_test_float("MC float");
DenseVectorSumTest<tags::CPU::MultiCore, double> mc_dense_vector_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorSumTest<tags::CPU::SSE, float> sse_dense_vector_sum_test_float("SSE float");
DenseVectorSumTest<tags::CPU::SSE, double> sse_dense_vector_sum_test_double("SSE double");
DenseVectorSumTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_sum_test_float("MC SSE float");
DenseVectorSumTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSumTest<tags::GPU::CUDA, float> cuda_dense_vector_sum_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSumTest<tags::Cell, float> cell_dense_vector_sum_test_float("Cell float");
DenseVectorSumTest<tags::Cell, double> cell_dense_vector_sum_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class DenseVectorSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sum_quick_test<" + type + ">")
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

            DenseVector<DataType_> dv1(size), dv2(size);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                *i = (i.index() + 1) / 1.23456789;
                *j = 2 * *i;
            }

            Sum<Tag_>::value(dv1, dv2);
            Sum<Tag_>::value(dv1, dv2);

            TEST(dv1.lock(lm_read_only),
                    for (typename Vector<DataType_>::ConstElementIterator i(dv1.begin_elements()),
                        i_end(dv1.end_elements()) ; i != i_end ; ++i)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_((i.index() + 1) * 5 / 1.23456789),
                        std::numeric_limits<DataType_>::epsilon() * 2 * (i.index() + 1) * 5 / 1.23456789);
                    },
                    dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv01(6, DataType_(1)), dv02(4, DataType_(1));
            TEST_CHECK_THROWS(Sum<Tag_>::value(dv01, dv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSumQuickTest<tags::CPU, float> dense_vector_sum_quick_test_float("float");
DenseVectorSumQuickTest<tags::CPU, double> dense_vector_sum_quick_test_double("double");
DenseVectorSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_sum_quick_test_float("MC float");
DenseVectorSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_sum_quick_test_float("SSE float");
DenseVectorSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_sum_quick_test_double("SSE double");
DenseVectorSumQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_sum_quick_test_float("MC SSE float");
DenseVectorSumQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSumQuickTest<tags::Cell, float> cell_dense_vector_sum_quick_test_float("Cell float");
DenseVectorSumQuickTest<tags::Cell, double> cell_dense_vector_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorSparseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSparseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sparse_vector_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
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

                Sum<Tag_>::value(dv1, sv2);

                TEST_CHECK_EQUAL(dv1, dv3);
            }

            DenseVector<DataType_> dv01(2);
            SparseVector<DataType_> sv02(1, 1);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dv01, sv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorSumTest<tags::CPU, float> dense_vector_sparse_vector_sum_test_float("float");
DenseVectorSparseVectorSumTest<tags::CPU, double> dense_vector_sparse_vector_sum_test_double("double");
DenseVectorSparseVectorSumTest<tags::CPU::MultiCore, float> mc_dense_vector_sparse_vector_sum_test_float("MC float");
DenseVectorSparseVectorSumTest<tags::CPU::MultiCore, double> mc_dense_vector_sparse_vector_sum_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorSparseVectorSumTest<tags::Cell, float> cell_dense_vector_sparse_vector_sum_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorSparseVectorSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSparseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sparse_vector_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
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

            Sum<Tag_>::value(dv1, sv2);

            TEST_CHECK_EQUAL(dv1, dv3);

            DenseVector<DataType_> dv01(2);
            SparseVector<DataType_> sv02(1, 1);
            TEST_CHECK_THROWS(Sum<Tag_>::value(dv01, sv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorSumQuickTest<tags::CPU, float> dense_vector_sparse_vector_sum_quick_test_float("float");
DenseVectorSparseVectorSumQuickTest<tags::CPU, double> dense_vector_sparse_vector_sum_quick_test_double("double");
DenseVectorSparseVectorSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_sparse_vector_sum_quick_test_float("MC float");
DenseVectorSparseVectorSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_sparse_vector_sum_quick_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorSparseVectorSumQuickTest<tags::Cell, float> cell_dense_vector_sparse_vector_sum_quick_test_float("Cell float");
#endif

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
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
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

