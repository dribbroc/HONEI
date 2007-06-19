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
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
 
#include <libla/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <string>
#include <tr1/memory>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class SparseVectorCopyTest :
    public BaseTest
{
    public:
        SparseVectorCopyTest(const std::string & type) :
            BaseTest("sparse_vector_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(20) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 10);
                std::tr1::shared_ptr<SparseVector<DataType_> > c(sv1.copy());

                for (typename Vector<DataType_>::ElementIterator i(c->begin_elements()), i_end(c->end_elements()) ;
                        i != i_end ; ++i)
                {
                    typename Vector<DataType_>::ConstElementIterator ci(i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ci, 0, std::numeric_limits<DataType_>::epsilon());

                    if (0 == (i.index() % 10))
                        *i = 1;
                }

                for (typename Vector<DataType_>::ConstElementIterator i(sv1.begin_elements()),
                        i_end(sv1.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
SparseVectorCopyTest<float> sparse_vector_copy_test_float("float");
SparseVectorCopyTest<double> sparse_vector_copy_test_double("double");

template <typename DataType_>
class SparseVectorCreationTest :
    public BaseTest
{
    public:
        SparseVectorCreationTest(const std::string & type) :
            BaseTest("sparse_vector_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > sv1(new SparseVector<DataType_>(size, size));
                std::tr1::shared_ptr<SparseVector<DataType_> > sv2(new SparseVector<DataType_>(size, 1));
                TEST_CHECK(true);
            }
        }
};
SparseVectorCreationTest<float> sparse_vector_creation_test_float("float");
SparseVectorCreationTest<double> sparse_vector_creation_test_double("double");
SparseVectorCreationTest<int> sparse_vector_creation_test_int("int");
SparseVectorCreationTest<bool> sparse_vector_creation_test_bool("bool");

template <typename DataType_>
class SparseVectorEqualityTest :
    public BaseTest
{
    public:
        SparseVectorEqualityTest(const std::string & type) :
            BaseTest("sparse_vector_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size);
                SparseVector<DataType_> sv2(size, size / 4 + 1);

                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()),
                        i_end(sv1.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                        *i = DataType_(1);
                }

                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()),
                        i_end(sv2.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                        *i = DataType_(1);
                }

                for (typename Vector<DataType_>::ConstElementIterator i(sv1.begin_elements()),
                        i_end(sv1.end_elements()), j(sv2.begin_elements()) ; i != i_end ;
                        ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }
                TEST_CHECK_EQUAL(sv1, sv2);
                
            }
        }
};
SparseVectorEqualityTest<float> sparse_vector_equality_test_float("float");
SparseVectorEqualityTest<double> sparse_vector_equality_test_double("double");

template <typename DataType_>
class SparseVectorFunctionsTest :
    public BaseTest
{
    public:
        SparseVectorFunctionsTest(const std::string & type) :
            BaseTest("sparse_vector_functions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size/3+1);

                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()),
                        i_end(sv1.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                        *i = DataType_(1);
                }
                
                TEST_CHECK_EQUAL(sv1.capacity(), size/3+1);
                TEST_CHECK_EQUAL(sv1.size(), size);
                TEST_CHECK_EQUAL(sv1.used_elements(), size/10);
            }
        }
};
SparseVectorFunctionsTest<float> sparse_vector_functions_test_float("float");
SparseVectorFunctionsTest<double> sparse_vector_functions_test_double("double");


template <typename DataType_>
class SparseVectorQuickTest :
    public QuickTest
{
    public:
        SparseVectorQuickTest(const std::string & type) :
            QuickTest("sparse_vector_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            SparseVector<DataType_> sv1(size, size);
            SparseVector<DataType_> sv2(size, 3);

            sv1[0]=static_cast<DataType_>(3.3598);
            sv2[0]=static_cast<DataType_>(3.3598);

            for (typename Vector<DataType_>::ConstElementIterator i(sv1.begin_elements()),
                    i_end(sv1.end_elements()), j(sv2.begin_elements()) ; i != i_end ;
                    ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
            }
            TEST_CHECK_EQUAL(sv1.capacity(), size);
            TEST_CHECK_EQUAL(sv2.capacity(), 3);
            TEST_CHECK_EQUAL(sv1.size(), size);
            TEST_CHECK_EQUAL(sv2.size(), size);            
            TEST_CHECK_EQUAL(sv1.used_elements(), 1);
            TEST_CHECK_EQUAL(sv2.used_elements(), 1); 
            TEST_CHECK_EQUAL(sv1, sv2);           
        }
};
SparseVectorQuickTest<float> sparse_vector_quick_test_float("float");
SparseVectorQuickTest<double> sparse_vector_quick_test_double("double");
