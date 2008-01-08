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
class ScalarBandedMatrixScaleTest :
    public BaseTest
{
    public:
        ScalarBandedMatrixScaleTest(const std::string & type) :
            BaseTest("scalar_banded_matrix_scale_test<" + type + ">")
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
ScalarBandedMatrixScaleTest<tags::CPU, float> scalar_banded_matrix_scale_test_float("float");
ScalarBandedMatrixScaleTest<tags::CPU, double> scalar_banded_matrix_scale_test_double("double");
ScalarBandedMatrixScaleTest<tags::CPU::MultiCore, float> mc_scalar_banded_matrix_scale_test_float("MC float");
ScalarBandedMatrixScaleTest<tags::CPU::MultiCore, double> mc_scalar_banded_matrix_scale_test_double("MC double");
#ifdef HONEI_SSE
ScalarBandedMatrixScaleTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_banded_matrix_scale_test_float("MC SSE float");
ScalarBandedMatrixScaleTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_banded_matrix_scale_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
//ScalarBandedMatrixScaleTest<tags::Cell, float> cell_scalar_banded_matrix_scale_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarBandedMatrixScaleQuickTest :
    public QuickTest
{
    public:
        ScalarBandedMatrixScaleQuickTest(const std::string & type) :
            QuickTest("scalar_banded_matrix_scale_quick_test<" + type + ">")
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
ScalarBandedMatrixScaleQuickTest<tags::CPU, float> scalar_banded_matrix_scale_quick_test_float("float");
ScalarBandedMatrixScaleQuickTest<tags::CPU, double> scalar_banded_matrix_scale_quick_test_double("double");
ScalarBandedMatrixScaleQuickTest<tags::CPU::MultiCore, float> mc_scalar_banded_matrix_scale_quick_test_float("MC float");
ScalarBandedMatrixScaleQuickTest<tags::CPU::MultiCore, double> mc_scalar_banded_matrix_scale_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarBandedMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_banded_matrix_scale_quick_test_float("MC SSE float");
ScalarBandedMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_banded_matrix_scale_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
//ScalarBandedMatrixScaleQuickTest<tags::Cell, float> cell_scalar_banded_matrix_scale_quick_test_float("Cell float");
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
DenseMatrixScaleTest<tags::CPU, float> scalar_dense_matrix_scale_test_float("float");
DenseMatrixScaleTest<tags::CPU, double> scalar_dense_matrix_scale_test_double("double");
DenseMatrixScaleTest<tags::CPU::MultiCore, float> mc_scalar_dense_matrix_scale_test_float("MC float");
DenseMatrixScaleTest<tags::CPU::MultiCore, double> mc_scalar_dense_matrix_scale_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixScaleTest<tags::CPU::SSE, float> sse_scalar_dense_matrix_scale_test_float("SSE float");
DenseMatrixScaleTest<tags::CPU::SSE, double> sse_scalar_dense_matrix_scale_test_double("SSE double");
DenseMatrixScaleTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_matrix_scale_test_float("MC SSE float");
DenseMatrixScaleTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_matrix_scale_test_double("MC SSE double");
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
class ScalarSparseMatrixScalarTest :
    public BaseTest
{
    public:
        ScalarSparseMatrixScalarTest(const std::string & type) :
            BaseTest("sparse_matrix_scale_test<" + type + ">")
        {
            register_tag(Tag_::name);
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
ScalarSparseMatrixScalarTest<tags::CPU, float> scalar_sparse_matrix_scalar_test_float("float");
ScalarSparseMatrixScalarTest<tags::CPU, double> scalar_sparse_matrix_scalar_test_double("double");
ScalarSparseMatrixScalarTest<tags::CPU::MultiCore, float> mc_scalar_sparse_matrix_scale_test_float("MC float");
ScalarSparseMatrixScalarTest<tags::CPU::MultiCore, double> mc_scalar_sparse_matrix_scale_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseMatrixScalarTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_matrix_scale_test_float("MC SSE float");
ScalarSparseMatrixScalarTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_matrix_scale_test_double("MC SSE double");
ScalarSparseMatrixScalarTest<tags::CPU::SSE, float> sse_scalar_sparse_matrix_scale_test_float("SSE float");
ScalarSparseMatrixScalarTest<tags::CPU::SSE, double> sse_scalar_sparse_matrix_scale_test_double("SSE double");
#endif
#ifdef HONEI_CELL
// NOT YET IMPLEMENTED:
//ScalarSparseMatrixScalarTest<tags::Cell, float> cell_scalar_sparse_matrix_scalar_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseMatrixScaleQuickTest :
    public QuickTest
{
    public:
        ScalarSparseMatrixScaleQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_scale_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
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
ScalarSparseMatrixScaleQuickTest<tags::CPU, float> scalar_sparse_matrix_scale_quick_test_float("float");
ScalarSparseMatrixScaleQuickTest<tags::CPU, double> scalar_sparse_matrix_scale_quick_test_double("double");
ScalarSparseMatrixScaleQuickTest<tags::CPU::MultiCore, float> mc_scalar_sparse_matrix_scalar_quick_test_float("MC float");
ScalarSparseMatrixScaleQuickTest<tags::CPU::MultiCore, double> mc_scalar_sparse_matrix_scalar_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_matrix_scalar_quick_test_float("MC SSE float");
ScalarSparseMatrixScaleQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_matrix_scalar_quick_test_double("MC SSE double");
ScalarSparseMatrixScaleQuickTest<tags::CPU::SSE, float> sse_scalar_sparse_matrix_scalar_quick_test_float("SSE float");
ScalarSparseMatrixScaleQuickTest<tags::CPU::SSE, double> sse_scalar_sparse_matrix_scalar_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
//ScalarSparseMatrixScalarQuickTest<tags::Cell, float> cell_scalar_sparse_matrix_scalar_quick_test_float("Cell float");
#endif



// Test cases for scaling vectors

template <typename Tag_, typename DataType_>
class ScalarDenseVectorScaleTest :
    public BaseTest
{
    public:
        ScalarDenseVectorScaleTest(const std::string & type) :
            BaseTest("dense_vector_scale_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> source(size * 4, DataType_(3));
                DenseVectorRange<DataType_> dv(source, size, 3);
                Scale<Tag_>::value(DataType_(2), dv);
                DataType_ v1(Norm<vnt_l_one>::value(dv));
                TEST_CHECK_EQUAL(v1, 6 * size);
            }
        }
};

ScalarDenseVectorScaleTest<tags::CPU, float> scalar_dense_vector_scale_test_float("float");
ScalarDenseVectorScaleTest<tags::CPU, double> scalar_dense_vector_scale_test_double("double");
ScalarDenseVectorScaleTest<tags::CPU::MultiCore, float> mc_scalar_dense_vector_scale_test_float("MC float");
ScalarDenseVectorScaleTest<tags::CPU::MultiCore, double> mc_scalar_dense_vector_scale_test_double("MC double");
#ifdef HONEI_SSE
ScalarDenseVectorScaleTest<tags::CPU::SSE, float> sse_scalar_dense_vector_scale_test_float("sse float");
ScalarDenseVectorScaleTest<tags::CPU::SSE, double> sse_scalar_dense_vector_scale_test_double("sse double");
ScalarDenseVectorScaleTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_vector_scale_test_float("MC SSE float");
ScalarDenseVectorScaleTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_vector_scale_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorScaleTest<tags::Cell, float> cell_scalar_dense_vector_scale_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarDenseVectorScaleQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorScaleQuickTest(const std::string & type) :
            QuickTest("dense_vector_scale_quick_test<" + type + ">")
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
            DenseVector<DataType_> source(size * 4, DataType_(3));
            DenseVectorRange<DataType_> dv1(source, size, 3);

            Scale<Tag_>::value(DataType_(2), dv1);
            DataType_ v1(Norm<vnt_l_one>::value(dv1));
            TEST_CHECK_EQUAL(v1, 6 * size);
        }
};
ScalarDenseVectorScaleQuickTest<tags::CPU, float>  scalar_dense_vector_scale_quick_test_float("float");
ScalarDenseVectorScaleQuickTest<tags::CPU, double> scalar_dense_vector_scale_quick_test_double("double");
ScalarDenseVectorScaleQuickTest<tags::CPU::MultiCore, float> mc_scalar_dense_vector_scale_quick_test_float("MC float");
ScalarDenseVectorScaleQuickTest<tags::CPU::MultiCore, double> mc_scalar_dense_vector_scale_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarDenseVectorScaleQuickTest<tags::CPU::SSE, float>  sse_scalar_dense_vector_scale_quick_test_float("sse float");
ScalarDenseVectorScaleQuickTest<tags::CPU::SSE, double> sse_scalar_dense_vector_scale_quick_test_double("sse double");
ScalarDenseVectorScaleQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_dense_vector_scale_quick_test_float("MC SSE float");
ScalarDenseVectorScaleQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_dense_vector_scale_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
ScalarDenseVectorScaleQuickTest<tags::Cell, float> cell_scalar_dense_vector_scale_quick_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseVectorScaleTest :
    public BaseTest
{
    public:
        ScalarSparseVectorScaleTest(const std::string & type) :
            BaseTest("sparse_vector_scale_test<" + type + ">")
        {
            register_tag(Tag_::name);
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

ScalarSparseVectorScaleTest<tags::CPU, float> scalar_sparse_vector_scale_test_float("float");
ScalarSparseVectorScaleTest<tags::CPU, double> scalar_sparse_vector_scale_test_double("double");
ScalarSparseVectorScaleTest<tags::CPU::MultiCore, float> mc_scalar_sparse_vector_scale_test_float("MC float");
ScalarSparseVectorScaleTest<tags::CPU::MultiCore, double> mc_scalar_sparse_vector_scale_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseVectorScaleTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_vector_scale_test_float("MC SSE float");
ScalarSparseVectorScaleTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_vector_scale_test_double("MC SSE double");
ScalarSparseVectorScaleTest<tags::CPU::SSE, float> sse_scalar_sparse_vector_scale_test_float("SSE float");
ScalarSparseVectorScaleTest<tags::CPU::SSE, double> sse_scalar_sparse_vector_scale_test_double("SSE double");
#endif
#ifdef HONEI_CELL
//ScalarSparseVectorScaleTest<tags::Cell, float> cell_scalar_sparse_vector_scale_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class ScalarSparseVectorScaleQuickTest :
    public QuickTest
{
    public:
        ScalarSparseVectorScaleQuickTest(const std::string & type) :
            QuickTest("sparse_vector_scale_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
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

ScalarSparseVectorScaleQuickTest<tags::CPU, float> scalar_sparse_vector_scale_quick_test_float("float");
ScalarSparseVectorScaleQuickTest<tags::CPU, double> scalar_sparse_vector_scale_quick_test_double("double");
ScalarSparseVectorScaleQuickTest<tags::CPU::MultiCore, float> mc_scalar_sparse_vector_scale_quick_test_float("MC float");
ScalarSparseVectorScaleQuickTest<tags::CPU::MultiCore, double> mc_scalar_sparse_vector_scale_quick_test_double("MC double");
#ifdef HONEI_SSE
ScalarSparseVectorScaleQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_scalar_sparse_vector_scale_quick_test_float("MC SSE float");
ScalarSparseVectorScaleQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_scalar_sparse_vector_scale_quick_test_double("MC SSE double");
ScalarSparseVectorScaleQuickTest<tags::CPU::SSE, float> sse_scalar_sparse_vector_scale_quick_test_float("SSE float");
ScalarSparseVectorScaleQuickTest<tags::CPU::SSE, double> sse_scalar_sparse_vector_scale_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
//ScalarSparseVectorScaleQuickTest<tags::Cell, float> cell_scalar_sparse_vector_scale_quick_test_float("Cell float");
#endif
