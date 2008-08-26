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

#include <honei/la/dense_vector_range.hh>
#include <honei/la/vector_error.hh>
#include <unittest/unittest.hh>

#include <cmath>
#include <string>


using namespace honei;
using namespace tests;

template <typename DataType_>
class DenseVectorRangeCreationTest :
    public BaseTest
{
    public:
        DenseVectorRangeCreationTest(const std::string & type) :
            BaseTest("dense_vector_range_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (2 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size, DataType_(0));
                DenseVectorRange<DataType_> dvr(dv, size - 1, 0);
                DenseVectorRange<DataType_> dvr_c(dvr, size - 3, 1);
                TEST_CHECK(true);
            }
        }
};
DenseVectorRangeCreationTest<float> dense_vector_range_creation_test_float("float");
DenseVectorRangeCreationTest<double> dense_vector_range_creation_test_double("double");


template <typename DataType_>
class DenseVectorRangeCopyTest :
    public BaseTest
{
    public:
        DenseVectorRangeCopyTest(const std::string & type) :
            BaseTest("dense_vector_range_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (2 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(1));
                DenseVectorRange<DataType_> dvr1(dv1, size - 1, 0), dvr2(dv2, size - 1, 0);
                std::tr1::shared_ptr<DenseVector<DataType_> > c(new DenseVector<DataType_>(dvr1.copy()));

                for (typename DenseVector<DataType_>::ElementIterator i(c->begin_elements()), i_end(c->end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename DenseVector<DataType_>::ConstElementIterator i(dvr1.begin_elements()),
                        i_end(dvr1.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }

                for (typename DenseVector<DataType_>::ConstElementIterator i(dvr2.begin_elements()),
                        i_end(dvr2.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 1, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseVectorRangeCopyTest<float> dense_vector_range_copy_test_float("float");
DenseVectorRangeCopyTest<double> dense_vector_range_copy_test_double("double");

template <typename DataType_>
class DenseVectorRangeEqualityTest :
    public BaseTest
{
    public:
        DenseVectorRangeEqualityTest(const std::string & type) :
            BaseTest("dense_vector_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (3 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv0(size, DataType_(1.23456));
                DenseVector<DataType_> dv1(size, DataType_(1.23456));
                DenseVectorRange<DataType_> dvr0(dv0, size - 1, 0);
                DenseVectorRange<DataType_> dvr1(dv1, size - 1, 1);

                for (typename DenseVectorRange<DataType_>::ElementIterator i(dvr0.begin_elements()), j(dvr1.begin_elements()),
                        i_end(dvr0.end_elements()) ; i != i_end ; ++i , ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }

                TEST_CHECK_EQUAL(dv0, dv1);

                DenseVectorRange<DataType_> dvr2(dv0, size, 0);
                DenseVectorRange<DataType_> dvr3(dvr0, size - 3, 2), dvr4(dvr1, size - 3, 1);
                for (typename DenseVectorRange<DataType_>::ElementIterator i(dvr3.begin_elements()), j(dvr4.begin_elements()),
                        i_end(dvr3.end_elements()) ; i != i_end; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }

                TEST_CHECK_EQUAL(dv0, dv1);

                TEST_CHECK_THROWS(dvr0 == dvr2, VectorSizeDoesNotMatch);
            }
        }
};
DenseVectorRangeEqualityTest<float> dense_vector_range_equality_test_float("float");
DenseVectorRangeEqualityTest<double> dense_vector_range_equality_test_double("double");

template <typename DataType_>
class DenseVectorRangeFunctionsTest :
    public BaseTest
{
    public:
        DenseVectorRangeFunctionsTest(const std::string & type) :
            BaseTest("dense_vector_range_functions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (2 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size);
                DenseVectorRange<DataType_> dvr(dv, size - 1, 0);
                for (typename DenseVector<DataType_>::ElementIterator i(dvr.begin_elements()), i_end(dvr.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_((i.index() + 1) / 1.23456789);
                }
                TEST_CHECK_EQUAL(dvr.size(), size - 1);

                for (int i = 0 ; i < size - 1 ; ++i)
                {
                    DataType_ s((i+1)/1.23456789);
                    TEST_CHECK_EQUAL_WITHIN_EPS((dvr)[i], s,
                            std::numeric_limits<DataType_>::epsilon());
                }
                for (int i = 0 ; i < size - 1 ; ++i)
                {
                    DataType_ s((i+5)/1.23456789);
                    (dvr)[i] = s;
                    TEST_CHECK_EQUAL_WITHIN_EPS((dvr)[i], s,
                        std::numeric_limits<DataType_>::epsilon());
                }

            }

            for (unsigned long size(10); size < (1 << 8); size <<=1)
            {
                DenseVector<DataType_> dv(size, DataType_(1.234));
                DenseVectorRange<DataType_> dvr(dv, size - 5, 3);

                for (typename DenseVectorRange<DataType_>::ElementIterator i(dvr.begin_elements()), i_end(dvr.end_elements());
                        i != i_end; ++i)
                {
                    *i = DataType_(3.98765421);
                }

                for (unsigned short i(0); i < 3; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS((dv)[i], DataType_(1.234), std::numeric_limits<DataType_>::epsilon());
                }
                for (typename DenseVector<DataType_>::ConstElementIterator ri(dv.element_at(3)), ri_end(dv.element_at(size - 2));
                        ri != ri_end; ++ri)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ri, DataType_(3.98765421), std::numeric_limits<DataType_>::epsilon());
                }
                for (unsigned long i(size - 2); i < size; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS((dv)[i], DataType_(1.234), std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseVectorRangeFunctionsTest<float> dense_vector_range_functions_test_float("float");
DenseVectorRangeFunctionsTest<double> dense_vector_range_functions_test_double("double");

template <typename DataType_>
class DenseVectorRangeQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseVector<DataType_> dv(4711, DataType_(123.987));
            DenseVectorRange<DataType_> dvr(dv, 238, 3101);
            DenseVectorRange<DataType_> dvr2(dvr, 113, 35);
            TEST_CHECK_EQUAL(dvr.size(), 238);
            TEST_CHECK_EQUAL(dvr, dvr);
            TEST_CHECK_EQUAL_WITHIN_EPS((dvr)[79] , 123.987, std::sqrt(std::numeric_limits<DataType_>::epsilon()));
            DataType_ s = DataType_(1.2345);
            (dvr)[39] = s;
            TEST_CHECK_EQUAL_WITHIN_EPS((dvr)[39] , s, sqrt(std::numeric_limits<DataType_>::epsilon()));
            TEST_CHECK_EQUAL_WITHIN_EPS((dvr2)[4] , s, sqrt(std::numeric_limits<DataType_>::epsilon()));

        #if 0
            DenseVectorRange<DataType_> dv2(10);
            for (int i = 0 ; i < 10 ; ++i)
            {
                dv2[i] = DataType_(i);
            }
            DenseVectorRange<DataType_> dv3(5);
            for (int i = 0 ; i < 5 ; ++i)
            {
                dv3[i] = DataType_(i + 3);
            }
            DenseVectorRange<DataType_> dv4(dv2, 3, 5);
            TEST_CHECK_EQUAL(dv3, dv4);
        #endif
        }
};
DenseVectorRangeQuickTest<float>  dense_vector_range_quick_test_float("float");
DenseVectorRangeQuickTest<double> dense_vector_range_quick_test_double("double");
