/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#include <libla/dense_vector_slice.hh>
#include <unittest/unittest.hh>

#include <string>


using namespace honei;
using namespace tests;

template <typename DataType_>
class DenseVectorSliceCreationTest :
    public BaseTest
{
    public:
        DenseVectorSliceCreationTest(const std::string & type) :
            BaseTest("dense_vector_slice_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size, DataType_(0));
                DenseVectorSlice<DataType_> dvs(dv, size - 5, 3, 1);
                TEST_CHECK(true);
            }
        }
};
DenseVectorSliceCreationTest<float> dense_vector_slice_creation_test_float("float");
DenseVectorSliceCreationTest<double> dense_vector_slice_creation_test_double("double");

template <typename DataType_>
class DenseVectorSliceCopyTest :
    public BaseTest
{
    public:
        DenseVectorSliceCopyTest(const std::string & type) :
            BaseTest("dense_vector_slice_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size, DataType_(0));
                DenseVectorSlice<DataType_> dvs(dv, size - 5, 3, 1);
                std::tr1::shared_ptr<DenseVector<DataType_> > c(new DenseVector<DataType_>(dvs.copy()));

                for (typename Vector<DataType_>::ElementIterator i(c->begin_elements()), i_end(c->end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename Vector<DataType_>::ConstElementIterator i(dvs.begin_elements()),
                        i_end(dvs.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseVectorSliceCopyTest<float> dense_vector_slice_copy_test_float("float");
DenseVectorSliceCopyTest<double> dense_vector_slice_copy_test_double("double");

template <typename DataType_>
class DenseVectorSliceEqualityTest :
    public BaseTest
{
    public:
        DenseVectorSliceEqualityTest(const std::string & type) :
            BaseTest("dense_vector_slice_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv0(size, DataType_(1.23456));

                for (typename Vector<DataType_>::ElementIterator i(dv0.begin_elements()), i_end(dv0.end_elements());
                        i != i_end; ++i)
                {
                    *i += i.index() % (size / 5);
                }

                DenseVectorSlice<DataType_> dvs0(dv0, size / 5, size / 5, 2);
                DenseVectorSlice<DataType_> dvs1(dv0, size / 5, size / 5 * 3, 2);

                for (typename Vector<DataType_>::ElementIterator i(dvs0.begin_elements()), j(dvs1.begin_elements()),
                    i_end(dvs0.end_elements()) ; i != i_end ; ++i , ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }

                TEST_CHECK_EQUAL(dvs0, dvs1);

                DenseVectorSlice<DataType_> dvs2(dv0, size / 2 - 2, 0, 2);
                TEST_CHECK_THROWS(dvs0 == dvs2, VectorSizeDoesNotMatch);
            }
        }
};
DenseVectorSliceEqualityTest<float> dense_vector_slice_equality_test_float("float");
DenseVectorSliceEqualityTest<double> dense_vector_slice_equality_test_double("double");

template <typename DataType_>
class DenseVectorSliceFunctionsTest :
    public BaseTest
{
    public:
        DenseVectorSliceFunctionsTest(const std::string & type) :
            BaseTest("dense_vector_slice_functions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size);
                DenseVectorSlice<DataType_> dvs(dv, size / 5, size / 5 * 2, 3);
                for (typename Vector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_((i.index() + 1) / 1.23456789);
                }
                TEST_CHECK_EQUAL(dvs.size(), size / 5);

                for (int i = 0 ; i < size / 5 ; ++i)
                {
                    DataType_ s((3 * i + size / 5 * 2 + 1) / 1.23456789);
                    TEST_CHECK_EQUAL_WITHIN_EPS(dvs[i], s,
                        std::numeric_limits<DataType_>::epsilon());
                }
                for (int i = 0 ; i < size / 5 ; ++i)
                {
                    DataType_ s((i + 5) / 1.23456789);
                    dvs[i] = s;
                    TEST_CHECK_EQUAL_WITHIN_EPS(dv[size / 5 * 2 + 3 * i], s,
                        std::numeric_limits<DataType_>::epsilon());
                }

            }
        }
};
DenseVectorSliceFunctionsTest<float> dense_vector_slice_functions_test_float("float");
DenseVectorSliceFunctionsTest<double> dense_vector_slice_functions_test_double("double");

template <typename DataType_>
class DenseVectorSliceQuickTest :
    public QuickTest
{
    public:
        DenseVectorSliceQuickTest(const std::string & type) :
            QuickTest("dense_vector_slice_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseVector<DataType_> dv(4711, DataType_(123.987));
            DenseVectorSlice<DataType_> dvs(dv, 238, 921, 13);
            TEST_CHECK_EQUAL(dv.size(), 4711);
            TEST_CHECK_EQUAL(dvs.size(), 238);
            TEST_CHECK_EQUAL(dv, dv);
            TEST_CHECK_EQUAL(dvs, dvs);
            TEST_CHECK_EQUAL_WITHIN_EPS(dv[4710] , 123.987, sqrt(std::numeric_limits<DataType_>::epsilon()));
            TEST_CHECK_EQUAL_WITHIN_EPS(dvs[89], 123.987, sqrt(std::numeric_limits<DataType_>::epsilon()));
            DataType_ s = DataType_(1.2345);
            dv[999] = s;
            dvs[67] = s + 5;
            TEST_CHECK_EQUAL_WITHIN_EPS(dv[999] , s, sqrt(std::numeric_limits<DataType_>::epsilon()));
            TEST_CHECK_EQUAL_WITHIN_EPS(dvs[67], s + 5, sqrt(std::numeric_limits<DataType_>::epsilon()));
            TEST_CHECK_EQUAL_WITHIN_EPS(dvs[6], s, sqrt(std::numeric_limits<DataType_>::epsilon()));
            TEST_CHECK_EQUAL_WITHIN_EPS(dv[1792], s + 5, sqrt(std::numeric_limits<DataType_>::epsilon()));
        }
};
DenseVectorSliceQuickTest<float>  dense_vector_slice_quick_test_float("float");
DenseVectorSliceQuickTest<double> dense_vector_slice_quick_test_double("double");
