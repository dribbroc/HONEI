/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#include <honei/la/dense_vector_range.hh>
#include <honei/la/norm.hh>
#include <honei/la/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;


template <typename Tag_, typename DataType_>
class DenseVectorNormValueTest :
    public BaseTest
{
    public:
        DenseVectorNormValueTest(const std::string & type) :
            BaseTest("dense_vector_norm_value_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_>dv(size);
                for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ s(size);

                DataType_ vmax(Norm<vnt_max, false, Tag_>::value(dv));
                DataType_ smax(s / 1.23456789);
                TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

                DataType_ v1(Norm<vnt_l_one>::value(dv));
                DataType_ s1(s * (s + 1) / 2 / 1.23456789);
                DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);

                DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(dv));
                DataType_ s2(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);

                DataType_ v3(Norm<vnt_l_two, true, Tag_>::value(dv));
                DataType_ s3(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                s3 = sqrt(s3);
                DataType_ eps3(s3 * 20 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v3, s3, eps3);
            }
        }
};

DenseVectorNormValueTest<tags::CPU, float> dense_vector_norm_value_test_float("float");
DenseVectorNormValueTest<tags::CPU, double> dense_vector_norm_value_test_double("double");
DenseVectorNormValueTest<tags::CPU::MultiCore, float> mc_dense_vector_norm_value_test_float("MC float");
DenseVectorNormValueTest<tags::CPU::MultiCore, double> mc_dense_vector_norm_value_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorNormValueTest<tags::Cell, float> cell_dense_vector_norm_value_test_float("float (Cell)");
DenseVectorNormValueTest<tags::Cell, double> cell_dense_vector_norm_value_test_double("double (Cell)");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorNormQuickTest :
    public QuickTest
{
    public:
        DenseVectorNormQuickTest(const std::string & type) :
            QuickTest("dense_vector_norm_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv(size);
            for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }

            DataType_ s(size);

            DataType_ vmax(Norm<vnt_max, false, Tag_>::value(dv));
            DataType_ smax(s / 1.23456789);
            TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

            DataType_ v1(Norm<vnt_l_one, false, Tag_>::value(dv));
            DataType_ s1(s * (s + 1) / 2 / 1.23456789);
            DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);

            DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(dv));
            DataType_ s2(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
            DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);

            DataType_ v3(Norm<vnt_l_two, true, Tag_>::value(dv));
            DataType_ s3(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
            s3 = sqrt(s3);
            DataType_ eps3(s3 * 20 * std::numeric_limits<float>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v3, s3, eps3);

        }
};
DenseVectorNormQuickTest<tags::CPU, float>  dense_vector_norm_quick_test_float("float");
DenseVectorNormQuickTest<tags::CPU, double> dense_vector_norm_quick_test_double("double");
DenseVectorNormQuickTest<tags::CPU::MultiCore, float>  mc_dense_vector_norm_quick_test_float("MC float");
DenseVectorNormQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_norm_quick_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorNormQuickTest<tags::Cell, float> cell_dense_vector_norm_quick_test_float("float (Cell)");
DenseVectorNormQuickTest<tags::Cell, double> cell_dense_vector_norm_quick_test_double("double (Cell)");
#endif


template <typename Tag_, typename DataType_>
class DenseVectorRangeNormValueTest :
    public BaseTest
{
    public:
        DenseVectorRangeNormValueTest(const std::string & type) :
            BaseTest("dense_vector_range_norm_value_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> d(size * 4, DataType_(100));
                for (int j(0) ; j < 4 ; j++)
                {
                    DenseVectorRange<DataType_>dv(d, size, j);
                    for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                            i != i_end ; ++i)
                    {
                        *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                    }

                    DataType_ s(size);

                    if ((Tag_::name != "sse") && (Tag_::name != "mc-sse"))
                    {
                        DataType_ vmax(Norm<vnt_max, false, Tag_>::value(dv));
                        DataType_ smax(s / 1.23456789);
                        TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

                        DataType_ v1(Norm<vnt_l_one>::value(dv));
                        DataType_ s1(s * (s + 1) / 2 / 1.23456789);
                        DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
                        TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
                    }

                    DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(dv));
                    DataType_ s2(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                    DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
                    TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);

                    DataType_ v3(Norm<vnt_l_two, true, Tag_>::value(dv));
                    DataType_ s3(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                    s3 = sqrt(s3);
                    DataType_ eps3(s3 * 20 * std::numeric_limits<float>::epsilon());
                    TEST_CHECK_EQUAL_WITHIN_EPS(v3, s3, eps3);
                }
            }
        }
};
DenseVectorRangeNormValueTest<tags::CPU, float> dense_vector_range_norm_value_test_float("float");
DenseVectorRangeNormValueTest<tags::CPU, double> dense_vector_range_norm_value_test_double("double");
DenseVectorRangeNormValueTest<tags::CPU::MultiCore, float> mc_dense_vector_range_norm_value_test_float("MC float");
DenseVectorRangeNormValueTest<tags::CPU::MultiCore, double> mc_dense_vector_range_norm_value_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeNormValueTest<tags::CPU::SSE, float> sse_dense_vector_range_norm_value_test_float("SSE float");
DenseVectorRangeNormValueTest<tags::CPU::SSE, double> sse_dense_vector_range_norm_value_test_double("SSE double");
DenseVectorRangeNormValueTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_range_norm_value_test_float("MC SSE float");
DenseVectorRangeNormValueTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_norm_value_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeNormValueTest<tags::GPU::CUDA, float> cuda_dense_vector_range_norm_value_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeNormValueTest<tags::Cell, float> cell_dense_vector_range_norm_value_test_float("float (Cell)");
DenseVectorRangeNormValueTest<tags::Cell, double> cell_dense_vector_range_norm_value_test_double("double (Cell)");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeNormQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeNormQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_norm_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> d(size * 4, DataType_(100));
            for (int j(0) ; j < 4 ; j++)
            {
                DenseVectorRange<DataType_> dv(d, size, j);
                for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ s(size);

                if ((Tag_::name != "sse") && (Tag_::name != "mc-sse"))
                {
                    DataType_ vmax(Norm<vnt_max, false, Tag_>::value(dv));
                    DataType_ smax(s / 1.23456789);
                    TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

                    DataType_ v1(Norm<vnt_l_one, false, Tag_>::value(dv));
                    DataType_ s1(s * (s + 1) / 2 / 1.23456789);
                    DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
                    TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
                }

                DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(dv));
                DataType_ s2(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);

                DataType_ v3(Norm<vnt_l_two, true, Tag_>::value(dv));
                DataType_ s3(s * (s + 1) * (2 * s + 1) / 6 / 1.23456789 / 1.23456789);
                s3 = sqrt(s3);
                DataType_ eps3(s3 * 20 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v3, s3, eps3);
            }
        }
};
DenseVectorRangeNormQuickTest<tags::CPU, float>  dense_vector_range_norm_quick_test_float("float");
DenseVectorRangeNormQuickTest<tags::CPU, double> dense_vector_range_norm_quick_test_double("double");
DenseVectorRangeNormQuickTest<tags::CPU::MultiCore, float>  mc_dense_vector_range_norm_quick_test_float("MC float");
DenseVectorRangeNormQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_norm_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeNormQuickTest<tags::CPU::SSE, float>  sse_dense_vector_range_norm_quick_test_float("SSE float");
DenseVectorRangeNormQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_norm_quick_test_double("SSE double");
DenseVectorRangeNormQuickTest<tags::CPU::MultiCore::SSE, float>  mc_sse_dense_vector_range_norm_quick_test_float("MC SSE float");
DenseVectorRangeNormQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_norm_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeNormQuickTest<tags::GPU::CUDA, float>  cuda_dense_vector_range_norm_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeNormQuickTest<tags::Cell, float> cell_dense_vector_range_norm_quick_test_float("float (Cell)");
DenseVectorRangeNormQuickTest<tags::Cell, double> cell_dense_vector_range_norm_quick_test_double("double (Cell)");
#endif


template <typename Tag_, typename DataType_>
class SparseVectorNormValueTest :
    public BaseTest
{
    public:
        SparseVectorNormValueTest(const std::string & type) :
            BaseTest("sparse_vector_norm_value_test<" + type + ">")
        {
            register_tag(Tag_::name); 
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv(size, size / 4 + 1);
                DataType_ smax(0);
                DataType_ s1(0);
                DataType_ s2(0);
                for (typename Vector<DataType_>::ElementIterator i(sv.begin_elements()),
                        i_end(sv.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_((i.index() + 1) / 1.23456789);
                        smax = DataType_((i.index() + 1) / 1.23456789);
                        s1 += DataType_((i.index() + 1) / 1.23456789);
                        s2 += DataType_((( i.index() + 1) / 1.23456789) * (( i.index() + 1) / 1.23456789));
                    }
                }

                DataType_ vmax(Norm<vnt_max, false, Tag_>::value(sv));
                TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

                DataType_ v1(Norm<vnt_l_one, false, Tag_>::value(sv));
                DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);

                DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(sv));
                DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);
            }
        }
};

SparseVectorNormValueTest<tags::CPU, float> sparse_vector_norm_value_test_float("float");
SparseVectorNormValueTest<tags::CPU, double> sparse_vector_norm_value_test_double("double");
SparseVectorNormValueTest<tags::CPU::MultiCore, float> mc_sparse_vector_norm_value_test_float("MC float");
SparseVectorNormValueTest<tags::CPU::MultiCore, double> mc_sparse_vector_norm_value_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseVectorNormQuickTest :
    public QuickTest
{
    public:
        SparseVectorNormQuickTest(const std::string & type) :
            QuickTest("sparse_vector_norm_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            SparseVector<DataType_> sv(size, size / 4 + 1);
            DataType_ smax(0);
            DataType_ s1(0);
            DataType_ s2(0);
            for (typename Vector<DataType_>::ElementIterator i(sv.begin_elements()),
                    i_end(sv.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_((i.index() + 1) / 1.23456789);
                    smax = DataType_((i.index() + 1) / 1.23456789);
                    s1 += DataType_((i.index() + 1) / 1.23456789);
                    s2 += DataType_((( i.index() + 1) / 1.23456789) * (( i.index() + 1) / 1.23456789));
                }
            }

            DataType_ vmax(Norm<vnt_max, false, Tag_>::value(sv));
            TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<float>::epsilon());

            DataType_ v1(Norm<vnt_l_one, false, Tag_>::value(sv));
            DataType_ eps1(s1 * 10 * std::numeric_limits<float>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);

            DataType_ v2(Norm<vnt_l_two, false, Tag_>::value(sv));
            DataType_ eps2(s2 * 20 * std::numeric_limits<float>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, eps2);
        }
};

SparseVectorNormQuickTest<tags::CPU, float> sparse_vector_norm_quick_test_float("float");
SparseVectorNormQuickTest<tags::CPU, double> sparse_vector_norm_quick_test_double("double");
SparseVectorNormQuickTest<tags::CPU::MultiCore, float> mc_sparse_vector_norm_quick_test_float("MC float");
SparseVectorNormQuickTest<tags::CPU::MultiCore, double> mc_sparse_vector_norm_quick_test_double("MC double");
