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
 
#include <libla/element_iterator.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class BandedMatrixElementIterationTest :
    public BaseTest
{
    public:
        BandedMatrixElementIterationTest(const std::string & type) :
            BaseTest("banded_matrix_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<BandedMatrix<DataType_> > bm(new BandedMatrix<DataType_>(size));

                typename Matrix<DataType_>::ConstElementIterator ce(bm->begin_elements()), ce_end(bm->end_elements());
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 0, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

BandedMatrixElementIterationTest<float> banded_matrix_element_iteration_test_float("float");
BandedMatrixElementIterationTest<double> banded_matrix_element_iteration_test_double("double");

template <typename DataType_>
class DenseMatrixElementIterationTest :
    public BaseTest
{
    public:
        DenseMatrixElementIterationTest(const std::string & type) :
            BaseTest("dense_matrix_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size,
                            static_cast<DataType_>(10)));

                typename MutableMatrix<DataType_>::ElementIterator e(dm->begin_elements()), e_end(dm->end_elements()) ;
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(e.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*e, 10, std::numeric_limits<DataType_>::epsilon());
                    *e = 333;
                    ++e;
                }

                for (typename MutableMatrix<DataType_>::ElementIterator me(dm->begin_elements()), me_end(dm->end_elements()) ;
                        me != me_end ; ++me)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*me, 333, std::numeric_limits<DataType_>::epsilon());
                }

                typename Matrix<DataType_>::ConstElementIterator ce(dm->begin_elements()), ce_end(dm->end_elements()) ;
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 333, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

DenseMatrixElementIterationTest<float> dense_matrix_element_iteration_test_float("float");
DenseMatrixElementIterationTest<double> dense_matrix_element_iteration_test_double("double");

template <typename DataType_>
class DenseVectorElementIterationTest :
    public BaseTest
{
    public:
        DenseVectorElementIterationTest(const std::string & type) :
            BaseTest("dense_vector_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size, static_cast<DataType_>(10)));

                typename Vector<DataType_>::ElementIterator e(dv->begin_elements()), e_end(dv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(e.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*e, 10, std::numeric_limits<DataType_>::epsilon());
                    *e = 222;
                    ++e;
                }

                typename Vector<DataType_>::ConstElementIterator ce(dv->begin_elements()), ce_end(dv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 222, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

DenseVectorElementIterationTest<float> dense_vector_element_iteration_test_float("float");
DenseVectorElementIterationTest<double> dense_vector_element_iteration_test_double("double");

template <typename DataType_>
class SparseVectorElementIterationTest :
    public BaseTest
{
    public:
        SparseVectorElementIterationTest(const std::string & type) :
            BaseTest("sparse_vector_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > 
                    sv(new SparseVector<DataType_>(size, (size / 5) + 1));
                
                typename Vector<DataType_>::ElementIterator f(sv->begin_elements()), f_end(sv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    if (i % 10 == 0) *f = 10;
                    ++f;
                }

                typename Vector<DataType_>::ElementIterator e(sv->begin_elements()), e_end(sv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(e.index(), i);
                    if (i % 10 == 0) 
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*e, 10, std::numeric_limits<DataType_>::epsilon());
                        *e = 222;
                    }
                    else TEST_CHECK_EQUAL_WITHIN_EPS(*e, 0, std::numeric_limits<DataType_>::epsilon());
                    ++e;
                }

                typename Vector<DataType_>::ConstElementIterator ce(sv->begin_elements()), ce_end(sv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    if (i % 10 == 0)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 222, std::numeric_limits<DataType_>::epsilon());
                    }
                    else TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 0, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
                
                unsigned long count(0);
                for (typename Vector<DataType_>::ConstElementIterator nz(sv->begin_non_zero_elements()),
                    nz_end(sv->end_non_zero_elements()) ; nz != nz_end ; ++nz )
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*nz, 222, std::numeric_limits<DataType_>::epsilon());
                    count ++;
                }
                TEST_CHECK_EQUAL(count, sv->used_elements());
            }
        }
};

SparseVectorElementIterationTest<float> sparse_vector_element_iteration_test_float("float");
SparseVectorElementIterationTest<double> sparse_vector_element_iteration_test_double("double");
