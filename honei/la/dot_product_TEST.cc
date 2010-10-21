/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
#include <honei/la/dot_product.hh>
#include <honei/la/norm.hh>
#include <honei/la/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <limits>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class DenseDotProductTest :
    public BaseTest
{
    public:
        DenseDotProductTest(const std::string & type) :
            BaseTest("dense_dot_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 13) ; size <<= 1)
            {
                DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

                DataType_ p0(DotProduct<Tag_>::value(dv1, dv0));
                DataType_ p1(DotProduct<Tag_>::value(dv1, dv1));
                TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());

                DenseVector<DataType_> dv2(size);
                for (typename DenseVector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v2(Norm<vnt_l_two, false>::value(dv2));
                DataType_ p2(DotProduct<Tag_>::value(dv2, dv2));
                //float eps(exp(-20 + 3.88127 * log(size)));
                //TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);
                TEST_CHECK_EQUAL_WITHIN_EPS((v2/p2), 1, (20*std::numeric_limits<DataType_>::epsilon()));
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};

DenseDotProductTest<tags::CPU, float> dense_scalar_product_test_float("float");
DenseDotProductTest<tags::CPU, double> dense_scalar_product_test_double("double");
DenseDotProductTest<tags::CPU::MultiCore, float> mc_dense_scalar_product_test_float("MC float");
DenseDotProductTest<tags::CPU::MultiCore, double> mc_dense_scalar_product_test_double("MC double");
#ifdef HONEI_SSE
DenseDotProductTest<tags::CPU::SSE, float> sse_dense_scalar_product_test_float("SSE float");
DenseDotProductTest<tags::CPU::SSE, double> sse_dense_scalar_product_test_double("SSE double");
DenseDotProductTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_scalar_product_test_float("MC SSE float");
DenseDotProductTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_scalar_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseDotProductTest<tags::GPU::CUDA, float> cuda_dense_scalar_product_test_float("float");
DenseDotProductTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_scalar_product_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseDotProductTest<tags::GPU::CUDA, double> cuda_dense_scalar_product_test_double("double");
DenseDotProductTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_scalar_product_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseDotProductQuickTest :
    public QuickTest
{
    public:
        DenseDotProductQuickTest(const std::string & type) :
            QuickTest("dense_dot_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(12345);
            DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

            DataType_ p0(DotProduct<Tag_>::value(dv1, dv0));
            DataType_ p1(DotProduct<Tag_>::value(dv1, dv1));
            TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());

            DenseVector<DataType_> dv2(size);
            for (typename DenseVector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_((i.index() + 1) / 1.23456789);
            }

            DataType_ v2(Norm<vnt_l_two, false>::value(dv2));
            DataType_ p2(DotProduct<Tag_>::value(dv2, dv2));
            //float eps(exp(-20 + 3.88127 * log(size))); // Cell limits...
            //TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS((v2/p2), 1, (5*std::numeric_limits<DataType_>::epsilon()));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseDotProductQuickTest<tags::CPU, float> dense_scalar_product_quick_test_float("float");
DenseDotProductQuickTest<tags::CPU, double> dense_scalar_product_quick_test_double("double");
DenseDotProductQuickTest<tags::CPU::MultiCore, float> mc_dense_scalar_product_quick_test_float("MC float");
DenseDotProductQuickTest<tags::CPU::MultiCore, double> mc_dense_scalar_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseDotProductQuickTest<tags::CPU::SSE, float> sse_dense_scalar_product_quick_test_float("SSE float");
DenseDotProductQuickTest<tags::CPU::SSE, double> sse_dense_scalar_product_quick_test_double("SSE double");
DenseDotProductQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_scalar_product_quick_test_float("MC SSE float");
DenseDotProductQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_scalar_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseDotProductQuickTest<tags::GPU::CUDA, float> cuda_dense_scalar_product_quick_test_float("float");
DenseDotProductQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_scalar_product_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseDotProductQuickTest<tags::GPU::CUDA, double> cuda_dense_scalar_product_quick_test_double("double");
DenseDotProductQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_scalar_product_quick_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseDotProductCellTest :
    public BaseTest
{
    public:
        DenseDotProductCellTest(const std::string & type) :
            BaseTest("dense_dot_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

                DataType_ p0(DotProduct<Tag_>::value(dv1, dv0));
                DataType_ p1(DotProduct<Tag_>::value(dv1, dv1));
                TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());

                DenseVector<DataType_> dv2(size);
                for (typename DenseVector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v2(Norm<vnt_l_two, false>::value(dv2));
                DataType_ p2(DotProduct<Tag_>::value(dv2, dv2));
                DataType_ eps(exp(-20 + 15 * log(log(size + 5))));
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};

#ifdef HONEI_CELL
DenseDotProductCellTest<tags::Cell, float> cell_dense_scalar_product_test_float("Cell float");
DenseDotProductCellTest<tags::Cell, double> cell_dense_scalar_product_test_double("Cell double");
#endif


template <typename Tag_, typename DataType_>
class DenseDotProductCellQuickTest :
    public QuickTest
{
    public:
        DenseDotProductCellQuickTest(const std::string & type) :
            QuickTest("dense_dot_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

            DataType_ p0(DotProduct<Tag_>::value(dv1, dv0));
            DataType_ p1(DotProduct<Tag_>::value(dv1, dv1));
            TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());

            DenseVector<DataType_> dv2(size);
            for (typename DenseVector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_((i.index() + 1) / 1.23456789);
            }

            DataType_ v2(Norm<vnt_l_two, false>::value(dv2));
            DataType_ p2(DotProduct<Tag_>::value(dv2, dv2));
            float eps(exp(-20 + 4 * log(size)));
            TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
#ifdef HONEI_CELL
DenseDotProductCellQuickTest<tags::Cell, float> cell_dense_dot_product_quick_test_float("Cell float");
DenseDotProductCellQuickTest<tags::Cell, double> cell_dense_dot_product_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeDotProductTest :
    public BaseTest
{
    public:
        DenseVectorRangeDotProductTest(const std::string & type) :
            BaseTest("dense_vector_range_dot_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 13) ; size <<= 1)
            {
                for (int i(0) ; i < 4 ; ++i)
                {
                    for (int j(0) ; j < 4 ; ++j)
                    {
                        DenseVector<DataType_> dv0 (size+3, DataType_(0)), dv1(size+3, DataType_(1));
                        DenseVectorRange<DataType_> dvr0 (dv0, size, i), dvr1(dv1, size, j);

                        DataType_ p0(DotProduct<Tag_>::value(dvr1, dvr0));
                        DataType_ p1(DotProduct<Tag_>::value(dvr1, dvr1));
                        TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                        TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());
                    }
                    DenseVector<DataType_> dv2(size+3);
                    for (typename DenseVector<DataType_>::ElementIterator k(dv2.begin_elements()), k_end(dv2.end_elements()) ;
                        k != k_end ; ++k)
                    {
                        *k = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                    }
                    DenseVectorRange<DataType_> dvr2 (dv2, size, i);
                    DataType_ v2(Norm<vnt_l_two, false>::value(dvr2));
                    DataType_ p2(DotProduct<Tag_>::value(dvr2, dvr2));
                    TEST_CHECK_EQUAL_WITHIN_EPS((v2/p2), 1, (20*std::numeric_limits<DataType_>::epsilon()));
                }
            }

            DenseVector<DataType_> dv00(4, DataType_(1)), dv01(4, DataType_(1));
            DenseVectorRange<DataType_> dvr00(dv00, 2, 2), dvr01(dv01, 3, 1);

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dvr00, dvr01), VectorSizeDoesNotMatch);
        }
};

DenseVectorRangeDotProductTest<tags::CPU, float> dense_vector_range_scalar_product_test_float("float");
DenseVectorRangeDotProductTest<tags::CPU, double> dense_vector_range_scalar_product_test_double("double");
DenseVectorRangeDotProductTest<tags::CPU::MultiCore, float> mc_dense_vector_range_scalar_product_test_float("MC float");
DenseVectorRangeDotProductTest<tags::CPU::MultiCore, double> mc_dense_vector_range_scalar_product_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeDotProductTest<tags::CPU::SSE, float> sse_dense_vector_range_scalar_product_test_float("SSE float");
DenseVectorRangeDotProductTest<tags::CPU::SSE, double> sse_dense_vector_range_scalar_product_test_double("SSE double");
DenseVectorRangeDotProductTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_scalar_product_test_float("MC SSE float");
DenseVectorRangeDotProductTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_scalar_product_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeDotProductTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scalar_product_test_float("float");
DenseVectorRangeDotProductTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_range_scalar_product_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeDotProductTest<tags::GPU::CUDA, double> cuda_dense_vector_range_scalar_product_test_double("double");
DenseVectorRangeDotProductTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_range_scalar_product_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeDotProductQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeDotProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_dot_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(12345);
            for (int i(0) ; i < 4 ; ++i)
            {
                for (int j(0) ; j < 4 ; ++j)
                {
                    DenseVector<DataType_> dv0 (size+3, DataType_(0)), dv1(size+3, DataType_(1));
                    DenseVectorRange<DataType_> dvr0 (dv0, size, i), dvr1(dv1, size, j);

                    DataType_ p0(DotProduct<Tag_>::value(dvr1, dvr0));
                    DataType_ p1(DotProduct<Tag_>::value(dvr1, dvr1));
                    TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                    TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());
                }
                DenseVector<DataType_> dv2(size+3);
                for (typename DenseVector<DataType_>::ElementIterator k(dv2.begin_elements()), k_end(dv2.end_elements()) ;
                    k != k_end ; ++k)
                {
                    *k = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                }
                DenseVectorRange<DataType_> dvr2 (dv2, size, i);
                DataType_ v2(Norm<vnt_l_two, false>::value(dvr2));
                DataType_ p2(DotProduct<Tag_>::value(dvr2, dvr2));
                TEST_CHECK_EQUAL_WITHIN_EPS((v2/p2), 1, (5*std::numeric_limits<DataType_>::epsilon()));
            }

            DenseVector<DataType_> dv00(4, DataType_(1)), dv01(4, DataType_(1));
            DenseVectorRange<DataType_> dvr00(dv00, 2, 2), dvr01(dv01, 3, 1);

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dvr00, dvr01), VectorSizeDoesNotMatch);
        }
};

DenseVectorRangeDotProductQuickTest<tags::CPU, float> dense_vector_range_scalar_product_quick_test_float("float");
DenseVectorRangeDotProductQuickTest<tags::CPU, double> dense_vector_range_scalar_product_quick_test_double("double");
DenseVectorRangeDotProductQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_range_scalar_product_quick_test_float("MC float");
DenseVectorRangeDotProductQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_scalar_product_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeDotProductQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_scalar_product_quick_test_float("SSE float");
DenseVectorRangeDotProductQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_scalar_product_quick_test_double("SSE double");
DenseVectorRangeDotProductQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_scalar_product_quick_test_float("MC SSE float");
DenseVectorRangeDotProductQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_scalar_product_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeDotProductQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scalar_product_quick_test_float("float");
DenseVectorRangeDotProductQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_range_scalar_product_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeDotProductQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_range_scalar_product_quick_test_double("double");
DenseVectorRangeDotProductQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_range_scalar_product_quick_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeDotProductCellTest :
    public BaseTest
{
    public:
        DenseVectorRangeDotProductCellTest(const std::string & type) :
            BaseTest("dense_vector_range_dot_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                for (int i(0) ; i < 4 ; ++i)
                {
                    for (int j(0) ; j < 4 ; ++j)
                    {
                        DenseVector<DataType_> dv0 (size+3, DataType_(0)), dv1(size+3, DataType_(1)), dv2(size+3, DataType_(1));
                        DenseVectorRange<DataType_> dvr0 (dv0, size, i), dvr1(dv1, size, j), dvr2(dv2, size, j);

                        DataType_ p0(DotProduct<Tag_>::value(dvr1, dvr0));
                        DataType_ p1(DotProduct<Tag_>::value(dvr1, dvr2));
                        TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                        TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());
                   }
                    DenseVector<DataType_> dv2(size+3), dv3(size+3);
                    for (typename DenseVector<DataType_>::ElementIterator k(dv2.begin_elements()), l(dv3.begin_elements()), k_end(dv2.end_elements()),
                            l_end(dv3.end_elements()) ; k != k_end ; ++k, ++l)
                    {
                        *k = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                        *l = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                    }
                    DenseVectorRange<DataType_> dvr2(dv2, size, i);
                    DenseVectorRange<DataType_> dvr3(dv3, size, i);

                    DataType_ v2(Norm<vnt_l_two, false>::value(dvr2));
                    DataType_ p2(DotProduct<Tag_>::value(dvr2, dvr3));
                    DataType_ eps(exp(-20 + 15 * log(log(size + 5))));
                    TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);
                }
            }

            DenseVector<DataType_> dv00(4, DataType_(1)), dv01(4, DataType_(1));
            DenseVectorRange<DataType_> dvr00(dv00, 2, 2), dvr01(dv01, 3, 1);

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dvr00, dvr01), VectorSizeDoesNotMatch);
        }
};

#ifdef HONEI_CELL
DenseVectorRangeDotProductCellTest<tags::Cell, float> cell_dense_vector_range_scalar_product_test_float("Cell float");
DenseVectorRangeDotProductCellTest<tags::Cell, double> cell_dense_vector_range_scalar_product_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeDotProductCellQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeDotProductCellQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_dot_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(120);
            for (int i(0) ; i < 4 ; ++i)
            {
                for (int j(0) ; j < 4 ; ++j)
                {
                    DenseVector<DataType_> dv0 (size+3, DataType_(0)), dv1(size+3, DataType_(1)), dv2(size+3, DataType_(1));
                    DenseVectorRange<DataType_> dvr0 (dv0, size, i), dvr1(dv1, size, j), dvr2(dv2, size, j);

                    DataType_ p0(DotProduct<Tag_>::value(dvr1, dvr0));
                    DataType_ p1(DotProduct<Tag_>::value(dvr1, dvr2));
                    TEST_CHECK_EQUAL_WITHIN_EPS(p0, 0, std::numeric_limits<DataType_>::epsilon());
                    TEST_CHECK_EQUAL_WITHIN_EPS(p1, size, std::numeric_limits<DataType_>::epsilon());
                }
                DenseVector<DataType_> dv2(size+3), dv3(size+3);
                for (typename DenseVector<DataType_>::ElementIterator k(dv2.begin_elements()), l(dv3.begin_elements()), k_end(dv2.end_elements()),
                        l_end(dv3.end_elements()) ; k != k_end ; ++k, ++l)
                {
                    *k = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                    *l = static_cast<DataType_>((k.index() + 1) / 1.23456789);
                }
                DenseVectorRange<DataType_> dvr2 (dv2, size, i);
                DenseVectorRange<DataType_> dvr3 (dv3, size, i);
                DataType_ v2(Norm<vnt_l_two, false>::value(dvr2));
                DataType_ p2(DotProduct<Tag_>::value(dvr2, dvr3));
                float eps(exp(-20 + 4 * log(size)));
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, eps);
            }

            DenseVector<DataType_> dv00(4, DataType_(1)), dv01(4, DataType_(1));
            DenseVectorRange<DataType_> dvr00(dv00, 2, 2), dvr01(dv01, 3, 1);

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(dvr00, dvr01), VectorSizeDoesNotMatch);
        }
};

#ifdef HONEI_CELL
DenseVectorRangeDotProductCellQuickTest<tags::Cell, float> cell_dense_vector_range_scalar_product_quick_test_float("Cell float");
DenseVectorRangeDotProductCellQuickTest<tags::Cell, double> cell_dense_vector_range_scalar_product_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class SparseDenseDotProductTest :
    public BaseTest
{
    public:
        SparseDenseDotProductTest(const std::string & type) :
            BaseTest("sparse_dense_dot_product_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DataType_ p1(0);
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                DenseVector<DataType_> dv2(size, DataType_(0));
                typename DenseVector<DataType_>::ElementIterator j(dv2.begin_elements());
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements());
                        i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        p1 += *i * *j;
                    }
                }

                DataType_ p0(DotProduct<Tag_>::value(sv1, dv2));
                DataType_ p2(DotProduct<Tag_>::value(dv2, sv1));
                TEST_CHECK_EQUAL_WITHIN_EPS(p0, p1, 2 * size * sqrt(std::numeric_limits<DataType_>::epsilon()));
                TEST_CHECK_EQUAL_WITHIN_EPS(p2, p1, 2 * size * sqrt(std::numeric_limits<DataType_>::epsilon()));
             }

            SparseVector<DataType_> sv00(1, 1);
            DenseVector<DataType_> dv01(2);

            TEST_CHECK_THROWS(DotProduct<>::value(sv00, dv01), VectorSizeDoesNotMatch);
            TEST_CHECK_THROWS(DotProduct<>::value(dv01, sv00), VectorSizeDoesNotMatch);
        }
};
SparseDenseDotProductTest<tags::CPU, float> sparse_dense_scalar_product_test_float("float");
SparseDenseDotProductTest<tags::CPU, double> sparse_dense_scalar_product_test_double("double");
//SparseDenseDotProductTest<tags::CPU::MultiCore, float> mc_sparse_dense_scalar_product_test_float("MC float");
//SparseDenseDotProductTest<tags::CPU::MultiCore, double> mc_sparse_dense_scalar_product_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseDenseDotProductQuickTest :
    public QuickTest
{
    public:
        SparseDenseDotProductQuickTest(const std::string & type) :
            QuickTest("sparse_dense_scalar_product_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DataType_ p1(0);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            DenseVector<DataType_> dv2(size, DataType_(0));
            typename DenseVector<DataType_>::ElementIterator j(dv2.begin_elements());
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                if (i.index() % 7 == 0)
                {
                    *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    p1 += *i * *j;
                }
            }

            DataType_ p0(DotProduct<>::value(sv1, dv2));
            TEST_CHECK_EQUAL(p0, p1);

            SparseVector<DataType_> sv00(1, 1);
            DenseVector<DataType_> dv01(2);

            TEST_CHECK_THROWS(DotProduct<Tag_>::value(sv00, dv01), VectorSizeDoesNotMatch);
        }
};
SparseDenseDotProductQuickTest<tags::CPU, float> sparse_dense_scalar_product_quick_test_float("float");
SparseDenseDotProductQuickTest<tags::CPU, double> sparse_dense_scalar_product_quick_test_double("double");
//SparseDenseDotProductQuickTest<tags::CPU::MultiCore, float> mc_sparse_dense_scalar_product_quick_test_float("MC float");
//SparseDenseDotProductQuickTest<tags::CPU::MultiCore, double> mc_sparse_dense_scalar_product_quick_test_double("MC double");

template <typename DataType_>
class SparseDotProductTest :
    public BaseTest
{
    public:
        SparseDotProductTest(const std::string & type) :
            BaseTest("sparse_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DataType_ p1(0);
                SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                        j(sv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                    }
                    if (i.index() % 7 == 0)
                    {
                        *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        p1 += *i * *j;
                    }
                }

                DataType_ p0(DotProduct<>::value(sv1, sv2));
                TEST_CHECK_EQUAL_WITHIN_EPS(p0, p1, std::numeric_limits<DataType_>::epsilon());
            }

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);

            TEST_CHECK_THROWS(DotProduct<>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseDotProductTest<float> sparse_scalar_product_test_float("float");
SparseDotProductTest<double> sparse_scalar_product_test_double("double");

template <typename DataType_>
class SparseDotProductQuickTest :
    public QuickTest
{
    public:
        SparseDotProductQuickTest(const std::string & type) :
            QuickTest("sparse_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            DataType_ p1(0);
            SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                if (i.index() % 7 == 0)
                {
                    *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    p1 += (*i) * (*j);
                }

                DataType_ p0(DotProduct<>::value(sv1, sv2));
                TEST_CHECK_EQUAL(p0, p1);
            }

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);

            TEST_CHECK_THROWS(DotProduct<>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseDotProductQuickTest<float> sparse_scalar_product_quick_test_float("float");
SparseDotProductQuickTest<double> sparse_scalar_product_quick_double("double");
