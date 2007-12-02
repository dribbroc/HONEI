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
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/norm.hh>
#include <libla/scale.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;

// Test cases for scaling matrices

template <typename Tag_, typename DataType_>
class ScalarBandedMatrixProductTest :
    public BaseTest
{
    public:
        ScalarBandedMatrixProductTest(const std::string & type) :
            BaseTest("scalar_banded_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                BandedMatrix<DataType_> bm1(size, dv1);
                Scale<Tag_>::value(DataType_(3), bm1);
                for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(bm1.begin_bands()),
                        ce_end(bm1.end_bands()) ; ce != ce_end ; ++ce)
                {
                    if (ce.index() != size - 1 )
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                            i_end(ce->end_elements()) ; i != i_end ; ++i)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                        }
                    }
                    else
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                            i_end(ce->end_elements()) ; i != i_end ; ++i)
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 6, std::numeric_limits<DataType_>::epsilon());
                        }
                    }
                }
            }
        }
};
ScalarBandedMatrixProductTest<tags::CPU, float> scalar_banded_matrix_product_test_float("float");
ScalarBandedMatrixProductTest<tags::CPU, double> scalar_banded_matrix_product_test_double("double");
ScalarBandedMatrixProductTest<tags::CPU::MultiCore, float> mc_scalar_banded_matrix_product_test_float("MC float");
ScalarBandedMatrixProductTest<tags::CPU::MultiCore, double> mc_scalar_banded_matrix_product_test_double("MC double");
#ifdef HONEI_SSE
ScalarBandedMatrixProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_banded_matrix_product_test_float("MC SSE float");
ScalarBandedMatrixProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_banded_matrix_product_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class ScalarBandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarBandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_banded_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> dv1(size, DataType_(2));
            BandedMatrix<DataType_> bm1(size, dv1);
            Scale<Tag_>::value(DataType_(3), bm1);
            for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(bm1.begin_bands()),
                    ce_end(bm1.end_bands()) ; ce != ce_end ; ++ce)
            {
                if (ce.index() != size - 1 )
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                            i_end(ce->end_elements()) ; i != i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    }
                }
                else
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(ce->begin_elements()),
                            i_end(ce->end_elements()) ; i != i_end ; ++i)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 6, std::numeric_limits<DataType_>::epsilon());
                    }
                }
            }
        }
};
ScalarBandedMatrixProductQuickTest<tags::CPU, float> scalar_banded_matrix_product_quick_test_float("float");
ScalarBandedMatrixProductQuickTest<tags::CPU, double> scalar_banded_matrix_product_quick_test_double("double");
ScalarBandedMatrixProductQuickTest<tags::CPU::MultiCore, float> mc_scalar_banded_matrix_product_quick_test_float("MC float");
ScalarBandedMatrixProductQuickTest<tags::CPU::MultiCore, double> mc_scalar_banded_matrix_product_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarBandedMatrixProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_banded_matrix_product_quick_test_float("MC SSE float");
ScalarBandedMatrixProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_banded_matrix_product_quick_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixScaleTest :
    public BaseTest
{
    public:
        DenseMatrixScaleTest(const std::string & type) :
            BaseTest("dense_matrix_scale_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(2)), dm2(size+1, size, DataType_(6));
                Scale<Tag_>::value(DataType_(3), dm1);

                TEST_CHECK_EQUAL(dm1, dm2);
            }
        }
};
DenseMatrixScaleTest<tags::CPU, float> scalar_dense_matrix_product_test_float("float");
DenseMatrixScaleTest<tags::CPU, double> scalar_dense_matrix_product_test_double("double");
DenseMatrixScaleTest<tags::CPU::MultiCore, float> mc_scalar_dense_matrix_product_test_float("MC float");
DenseMatrixScaleTest<tags::CPU::MultiCore, double> mc_scalar_dense_matrix_product_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixScaleTest<tags::CPU::SSE, float> sse_scalar_dense_matrix_product_test_float("SSE float");
DenseMatrixScaleTest<tags::CPU::SSE, double> sse_scalar_dense_matrix_product_test_double("SSE double");
DenseMatrixScaleTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_matrix_product_test_float("MC SSE float");
DenseMatrixScaleTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_matrix_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixScaleTest<tags::Cell, float> cell_dense_matrix_scale_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixScaleQuickTest :
    public QuickTest
{
    public:
        DenseMatrixScaleQuickTest(const std::string & type) :
            QuickTest("dense_matrix_scale_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(100);
            DenseMatrix<DataType_> dm(size+1, size, DataType_(2));
            Scale<Tag_>::value(DataType_(3), dm);

            DataType_ vsum(0), ssum(2 * DataType_(size) * DataType_(size + 1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                vsum += *i;
            }
            TEST_CHECK_EQUAL_WITHIN_EPS(vsum, static_cast<DataType_>(6 * size * (size +1)),
                std::numeric_limits<DataType_>::epsilon());
        }
};
DenseMatrixScaleQuickTest<tags::CPU, float> dense_matrix_scale_quick_test_float("float");
DenseMatrixScaleQuickTest<tags::CPU, double> dense_matrix_scale_quick_test_double("double");
DenseMatrixScaleQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_scale_quick_test_float("MC float");
DenseMatrixScaleQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_scale_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixScaleQuickTest<tags::CPU::SSE, float> sse_dense_matrix_scale_quick_test_float("SSE float");
DenseMatrixScaleQuickTest<tags::CPU::SSE, double> sse_dense_matrix_scale_quick_test_double("SSE double");
DenseMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_matrix_scale_quick_test_float("MC SSE float");
DenseMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_matrix_scale_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixScaleQuickTest<tags::Cell, float> cell_dense_matrix_scale_quick_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseMatrixProductTest :
    public BaseTest
{
    public:
        ScalarSparseMatrixProductTest(const std::string & type) :
            BaseTest("scalar_sparse_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                    sm2(size+1, size, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(2);
                        *j = DataType_(6);
                    }
                }
                Scale<Tag_>::value(DataType_(3), sm1);

                TEST_CHECK_EQUAL(sm1, sm2);
            }
        }
};
ScalarSparseMatrixProductTest<tags::CPU, float> scalar_sparse_matrix_product_test_float("float");
ScalarSparseMatrixProductTest<tags::CPU, double> scalar_sparse_matrix_product_test_double("double");
ScalarSparseMatrixProductTest<tags::CPU::MultiCore, float> mc_scalar_sparse_matrix_product_test_float("MC float");
ScalarSparseMatrixProductTest<tags::CPU::MultiCore, double> mc_scalar_sparse_matrix_product_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseMatrixProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_matrix_product_test_float("MC SSE float");
ScalarSparseMatrixProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_matrix_product_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarSparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_sparse_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                sm2(size+1, size, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(2);
                    *j = DataType_(6);
                }
            }
            Scale<Tag_>::value(DataType_(3), sm1);

            TEST_CHECK_EQUAL(sm1, sm2);
        }
};
ScalarSparseMatrixProductQuickTest<tags::CPU, float> scalar_sparse_matrix_product_quick_test_float("float");
ScalarSparseMatrixProductQuickTest<tags::CPU, double> scalar_sparse_matrix_product_quick_test_double("double");
ScalarSparseMatrixProductQuickTest<tags::CPU::MultiCore, float> mc_scalar_sparse_matrix_product_quick_test_float("MC float");
ScalarSparseMatrixProductQuickTest<tags::CPU::MultiCore, double> mc_scalar_sparse_matrix_product_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_matrix_product_quick_test_float("MC SSE float");
ScalarSparseMatrixProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_matrix_product_quick_test_double("MC SSE double");
#endif

// Test cases for scaling vectors

template <typename Tag_, typename DataType_>
class ScalarDenseVectorProductTest :
    public BaseTest
{
    public:
        ScalarDenseVectorProductTest(const std::string & type) :
            BaseTest("dense_vector_scale_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));

                Scale<Tag_>::value(DataType_(2), dv1);
                DataType_ v1(Norm<vnt_l_one>::value(dv1));
                TEST_CHECK_EQUAL(v1, 6 * size);
            }
        }
};

ScalarDenseVectorProductTest<tags::CPU, float> scalar_dense_vector_product_test_float("float");
ScalarDenseVectorProductTest<tags::CPU, double> scalar_dense_vector_product_test_double("double");
ScalarDenseVectorProductTest<tags::CPU::MultiCore, float> mc_scalar_dense_vector_product_test_float("MC float");
ScalarDenseVectorProductTest<tags::CPU::MultiCore, double> mc_scalar_dense_vector_product_test_double("MC double");
#ifdef HONEI_SSE
ScalarDenseVectorProductTest<tags::CPU::SSE, float> sse_scalar_dense_vector_product_test_float("sse float");
ScalarDenseVectorProductTest<tags::CPU::SSE, double> sse_scalar_dense_vector_product_test_double("sse double");
ScalarDenseVectorProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_vector_product_test_float("MC SSE float");
ScalarDenseVectorProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_vector_product_test_double("MC SSE double");
#endif

template <typename Tag_, typename DataType_>
class ScalarDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_scale_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(3));

            Scale<Tag_>::value(DataType_(2), dv1);
            DataType_ v1(Norm<vnt_l_one>::value(dv1));
            TEST_CHECK_EQUAL(v1, 6 * size);
        }
};
ScalarDenseVectorProductQuickTest<tags::CPU, float>  scalar_dense_vector_product_quick_test_float("float");
ScalarDenseVectorProductQuickTest<tags::CPU, double> scalar_dense_vector_product_quick_test_double("double");
ScalarDenseVectorProductQuickTest<tags::CPU::MultiCore, float> mc_scalar_dense_vector_product_quick_test_float("MC float");
ScalarDenseVectorProductQuickTest<tags::CPU::MultiCore, double> mc_scalar_dense_vector_product_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarDenseVectorProductQuickTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_product_quick_test_float("sse float");
ScalarDenseVectorProductQuickTest<tags::CPU::SSE, double> sse_scalar_dense_vector_product_quick_test_double("sse double");
ScalarDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_vector_product_quick_test_float("MC SSE float");
ScalarDenseVectorProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_vector_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorProductQuickTest<tags::Cell, float> cell_scalar_dense_vector_product_quick_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseVectorProductTest :
    public BaseTest
{
    public:
        ScalarSparseVectorProductTest(const std::string & type) :
            BaseTest("scalar_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = 3;
                }
                Scale<Tag_>::value(static_cast<DataType_>(2), sv1);
                DataType_ v1(Norm<vnt_l_one>::value(sv1));
                TEST_CHECK_EQUAL(v1, 6 * (size / 10 + 1));
            }
        }
};

ScalarSparseVectorProductTest<tags::CPU, float> scalar_sparse_vector_product_test_float("float");
ScalarSparseVectorProductTest<tags::CPU, double> scalar_sparse_vector_product_test_double("double");
ScalarSparseVectorProductTest<tags::CPU::MultiCore, float> mc_scalar_sparse_vector_product_test_float("MC float");
ScalarSparseVectorProductTest<tags::CPU::MultiCore, double> mc_scalar_sparse_vector_product_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseVectorProductTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_vector_product_test_float("MC SSE float");
ScalarSparseVectorProductTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_vector_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
ScalarSparseVectorProductTest<tags::Cell, float> cell_scalar_sparse_vector_product_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        ScalarSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("scalar_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(111);
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = 3;
            }
            Scale<Tag_>::value(static_cast<DataType_>(2), sv1);
            DataType_ v1(Norm<vnt_l_one>::value(sv1));
            TEST_CHECK_EQUAL(v1, 6 * (size / 10 + 1));
        }
};

ScalarSparseVectorProductQuickTest<tags::CPU, float> scalar_sparse_vector_product_quick_test_float("float");
ScalarSparseVectorProductQuickTest<tags::CPU, double> scalar_sparse_vector_product_quick_test_double("double");
ScalarSparseVectorProductQuickTest<tags::CPU::MultiCore, float> mc_scalar_sparse_vector_product_quick_test_float("MC float");
ScalarSparseVectorProductQuickTest<tags::CPU::MultiCore, double> mc_scalar_sparse_vector_product_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseVectorProductQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_vector_product_quick_test_float("MC SSE float");
ScalarSparseVectorProductQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_vector_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
ScalarSparseVectorProductQuickTest<tags::Cell, float> cell_scalar_sparse_vector_product_quick_test_float("Cell float");
#endif
