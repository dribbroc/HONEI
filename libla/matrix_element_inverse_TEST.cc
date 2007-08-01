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

#include <libla/matrix_element_inverse.hh>
#include <unittest/unittest.hh>

#include <iostream>
#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class BandedMatrixElementInverseTest :
    public BaseTest
{
    public:
        BandedMatrixElementInverseTest(const std::string & type) :
            BaseTest("banded_matrix_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size));                    
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size)); 
                for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()),
                        j(dv2->begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);  
                
                TEST_CHECK_EQUAL(MatrixElementInverse<>::value(bm1), bm2);
            }
        }
};
BandedMatrixElementInverseTest<float> banded_matrix_element_inverse_test_float("float");
BandedMatrixElementInverseTest<double> banded_matrix_element_inverse_test_double("double");

template <typename DataType_>
class BandedMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("banded_matrix_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size));                    
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size)); 
                for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()),
                        j(dv2->begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);  
                
                TEST_CHECK_EQUAL(MatrixElementInverse<>::value(bm1), bm2);
            }
        }
};
BandedMatrixElementInverseQuickTest<float> banded_matrix_element_inverse_quick_test_float("float");
BandedMatrixElementInverseQuickTest<double> banded_matrix_element_inverse_quick_test_double("double");

template <typename DataType_>
class DenseMatrixElementInverseTest :
    public BaseTest
{
    public:
        DenseMatrixElementInverseTest(const std::string & type) :
            BaseTest("dense_matrix_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(size, size,
                            static_cast<DataType_>(0))),
                        dm2(new DenseMatrix<DataType_>(size, size, static_cast<DataType_>(0)));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1->begin_elements()), i_end(dm1->end_elements()),
                        j(dm2->begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }

                TEST_CHECK_EQUAL(MatrixElementInverse<>::value(*dm1), *dm2);
            }
        }
};
DenseMatrixElementInverseTest<float> dense_matrix_element_inverse_test_float("float");
DenseMatrixElementInverseTest<double> dense_matrix_element_inverse_test_double("double");

template <typename DataType_>
class DenseMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_matrix_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(3, 2,
                        static_cast<DataType_>(0))),
                    dm2(new DenseMatrix<DataType_>(3, 2, static_cast<DataType_>(0)));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1->begin_elements()), i_end(dm1->end_elements()),
                    j(dm2->begin_elements()) ; i != i_end ; ++i)
            {
                *i = i.index() + 1;
                *j = 1 / static_cast<DataType_>(i.index() + 1);
                ++j;
            }

            TEST_CHECK_EQUAL(MatrixElementInverse<>::value(*dm1), *dm2);
        }
};
DenseMatrixElementInverseQuickTest<float>  dense_matrix_element_inverse_quick_test_float("float");
DenseMatrixElementInverseQuickTest<double> dense_matrix_element_inverse_quick_test_double("double");

template <typename DataType_>
class SparseMatrixElementInverseTest :
    public BaseTest
{
    public:
        SparseMatrixElementInverseTest(const std::string & type) :
            BaseTest("sparse_matrix_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                    sm2(size, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()), j(sm2.begin_elements()); 
                    i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = i.index() + 1;
                        *j = 1 / static_cast<DataType_>(i.index() + 1);

                    }
                }  
                TEST_CHECK_EQUAL(MatrixElementInverse<>::value(sm1), sm2);
            }
        }
};
SparseMatrixElementInverseTest<float> sparse_matrix_element_inverse_test_float("float");
SparseMatrixElementInverseTest<double> sparse_matrix_element_inverse_test_double("double");

template <typename DataType_>
class SparseMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()), j(sm2.begin_elements()); 
                i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);

                }
            }  
            TEST_CHECK_EQUAL(MatrixElementInverse<>::value(sm1), sm2);
        }
};
SparseMatrixElementInverseQuickTest<float> sparse_matrix_element_inverse_quick_test_float("float");
SparseMatrixElementInverseQuickTest<double> sparse_matrix_element_inverse_quick_test_double("double");
