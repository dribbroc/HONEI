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

#include <libla/vector_element_inverse.hh>
#include <unittest/unittest.hh>

#include <iostream>
#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorElementInverseTest :
    public BaseTest
{
    public:
        DenseVectorElementInverseTest(const std::string & type) :
            BaseTest("vector_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > 
                    dv1(new DenseVector<DataType_>(size, static_cast<DataType_>(0))),
                    dv2(new DenseVector<DataType_>(size, static_cast<DataType_>(0)));
                for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()),
                        j(dv2->begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }

                TEST_CHECK_EQUAL(VectorElementInverse<>::value(*dv1), *dv2);
            }
        }
};

DenseVectorElementInverseTest<float> dense_vector_element_inverse_test_float("float");
DenseVectorElementInverseTest<double> dense_vector_element_inverse_test_double("double");

template <typename DataType_>
class DenseVectorElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_vector_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<DenseVector<DataType_> > 
                dv1(new DenseVector<DataType_>(3, static_cast<DataType_>(0))),
                dv2(new DenseVector<DataType_>(3, static_cast<DataType_>(0)));
            for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()),
                    j(dv2->begin_elements()) ; i != i_end ; ++i)
            {
                *i = i.index() + 1;
                *j = 1 / static_cast<DataType_>(i.index() + 1);
                ++j;
            }

            TEST_CHECK_EQUAL(VectorElementInverse<>::value(*dv1), *dv2);
        }
};
DenseVectorElementInverseQuickTest<float>  dense_vector_element_inverse_quick_test_float("float");
DenseVectorElementInverseQuickTest<double> dense_vector_element_inverse_quick_test_double("double");


template <typename DataType_>
class SparseVectorElementInverseTest :
    public BaseTest
{
    public:
        SparseVectorElementInverseTest(const std::string & type) :
            BaseTest("sparse_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > 
                    sv1(new SparseVector<DataType_>(size, size / 7 + 1)),
                    sv2(new SparseVector<DataType_>(size, size / 8 + 1));
                for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()),
                        j(sv2->begin_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = i.index() + 1;
                        *j = 1 / static_cast<DataType_>(i.index() + 1);
                    }
                    ++j;
                }

                TEST_CHECK_EQUAL(VectorElementInverse<>::value(*sv1), *sv2);
            }
        }
};
SparseVectorElementInverseTest<float> sparse_vector_element_inverse_test_float("float");
SparseVectorElementInverseTest<double> sparse_vector_element_inverse_test_double("double");

template <typename DataType_>
class SparseVectorElementInverseQuickTest :
    public QuickTest
{
    public:
        SparseVectorElementInverseQuickTest(const std::string & type) :
            QuickTest("sparse_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (20);
            std::tr1::shared_ptr<SparseVector<DataType_> > 
                sv1(new SparseVector<DataType_>(size, size / 7 + 1)),
                sv2(new SparseVector<DataType_>(size, size / 8 + 1));
            for (typename Vector<DataType_>::ElementIterator i(sv1->begin_elements()), i_end(sv1->end_elements()),
                    j(sv2->begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                }
                ++j;
            }

            TEST_CHECK_EQUAL(VectorElementInverse<>::value(*sv1), *sv2);
        }
};
SparseVectorElementInverseQuickTest<float> sparse_vector_element_inverse_quick_test_float("float");
SparseVectorElementInverseQuickTest<double> sparse_vector_element_inverse_quick_test_double("double");
