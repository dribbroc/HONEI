/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
 * Copyright (c) 2009 Sven Mallach <mallach@honei.org>
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

#include <honei/la/dense_vector_slice.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/unittest.hh>

#include <cmath>
#include <cstdlib>
#include <limits>
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
                TEST_CHECK_EQUAL(*(dvs.begin_elements()), *(dv.element_at(3)));
                TEST_CHECK_EQUAL(*(dvs.end_elements()), *(dv.element_at(3 + size - 5)));
            }
        }
};
DenseVectorSliceCreationTest<float> dense_vector_slice_creation_test_float("float");
DenseVectorSliceCreationTest<double> dense_vector_slice_creation_test_double("double");

template <typename DataType_>
class DenseVectorSliceIterationTest :
    public BaseTest
{
    public:
        DenseVectorSliceIterationTest(const std::string & type) :
            BaseTest("dense_vector_slice_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size);
                for (typename DenseVector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index();
                }

                unsigned pos(rand() % size / 2);

                DenseVectorSlice<DataType_> dvs(dv, size / 4, pos, 2);

                typename DenseVector<DataType_>::ConstElementIterator j(dv.element_at(pos));
                for (typename DenseVectorSlice<DataType_>::ConstElementIterator i(dvs.begin_elements()), i_end(dvs.end_elements()) ;
                        i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL(*i, *j);
                    ++j;
                }
            }
        }
};
DenseVectorSliceIterationTest<float> dense_vector_slice_iteration_test_float("float");
DenseVectorSliceIterationTest<double> dense_vector_slice_itearation_test_double("double");

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
                shared_ptr<DenseVector<DataType_> > c(new DenseVector<DataType_>(dvs.copy()));

                for (typename DenseVector<DataType_>::ElementIterator i(c->begin_elements()), i_end(c->end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename DenseVectorSlice<DataType_>::ConstElementIterator i(dvs.begin_elements()),
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

                for (typename DenseVector<DataType_>::ElementIterator i(dv0.begin_elements()), i_end(dv0.end_elements());
                        i != i_end; ++i)
                {
                    *i += i.index() % (size / 5);
                }

                DenseVectorSlice<DataType_> dvs0(dv0, size / 5, size / 5, 2);
                DenseVectorSlice<DataType_> dvs1(dv0, size / 5, size / 5 * 3, 2);

                for (typename DenseVectorSlice<DataType_>::ElementIterator i(dvs0.begin_elements()), j(dvs1.begin_elements()),
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
                for (typename DenseVector<DataType_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_((i.index() + 1) / 1.23456789);
                }
                TEST_CHECK_EQUAL(dvs.size(), size / 5);

                for (unsigned long i = 0 ; i < size / 5 ; ++i)
                {
                    DataType_ s((3 * i + size / 5 * 2 + 1) / 1.23456789);
                    TEST_CHECK_EQUAL_WITHIN_EPS(dvs[i], s,
                        std::numeric_limits<DataType_>::epsilon());
                }
                for (unsigned long i = 0 ; i < size / 5 ; ++i)
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
