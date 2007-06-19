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
 
#include <libla/dense_matrix.hh>
#include <libla/matrix_mask.hh>
#include <libla/matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseMatrixMaskTest :
    public BaseTest
{
    public:
        DenseMatrixMaskTest(const std::string & type) :
            BaseTest("dense_matrix_mask_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                std::tr1::shared_ptr<DenseMatrix<bool> > mask(new DenseMatrix<bool>(
                    size, size+1, false));
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(
                    size, size+1, static_cast<DataType_>(1)));
                
                for (typename MutableMatrix<bool>::ElementIterator i(mask->begin_elements()), 
                    i_end(mask->end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index()%2 == 0) 
                    {
                        *i=true;
                        ++count;
                    }
                    
                }
                DenseMatrix<DataType_> result(MatrixMask<DataType_>::value(*dm, *mask));
                DataType_ sum = 0;
                for (typename MutableMatrix<DataType_>::ElementIterator j(result.begin_elements()), 
                    j_end(result.end_elements()) ; j != j_end ; ++j)
                {
                    sum += *j;
                }
                
                TEST_CHECK_EQUAL(sum, count);
            }

            std::tr1::shared_ptr<Matrix<bool> > mask1(new DenseMatrix<bool>(
                2, 3, false));
            std::tr1::shared_ptr<Matrix<bool> > mask2(new DenseMatrix<bool>(
                3, 4, false));                
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(
                2, 4, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask1), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask2), MatrixColumnsDoNotMatch);
        }
};

DenseMatrixMaskTest<float> dense_matrix_mask_test_float("float");
DenseMatrixMaskTest<double> dense_matrix_mask_test_double("double"); 

template <typename DataType_>
class DenseMatrixMaskQuickTest :
    public QuickTest
{
    public:
        DenseMatrixMaskQuickTest(const std::string & type) :
            QuickTest("dense_matrix_mask_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DataType_ count = 0;
            std::tr1::shared_ptr<DenseMatrix<bool> > mask(new DenseMatrix<bool>(
                size, size+1, false));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(
                size, size+1, static_cast<DataType_>(1)));
            
            for (typename MutableMatrix<bool>::ElementIterator i(mask->begin_elements()), 
                i_end(mask->end_elements()) ; i != i_end ; ++i)
            {
                if (i.index()%2 == 0) 
                {
                    *i=true;
                    ++count;
                }
                
            }
            DenseMatrix<DataType_> result(MatrixMask<DataType_>::value(*dm, *mask));
            DataType_ sum = 0;
            for (typename MutableMatrix<DataType_>::ElementIterator j(result.begin_elements()), 
                j_end(result.end_elements()) ; j != j_end ; ++j)
            {
                sum += *j;
            }
            
            TEST_CHECK_EQUAL(sum, count);

            std::tr1::shared_ptr<Matrix<bool> > mask1(new DenseMatrix<bool>(
                2, 3, false));
            std::tr1::shared_ptr<Matrix<bool> > mask2(new DenseMatrix<bool>(
                3, 4, false));                
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(
                2, 4, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask1), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask2), MatrixColumnsDoNotMatch);

        }
};
DenseMatrixMaskQuickTest<float>  dense_matrix_mask_quick_test_float("float");
DenseMatrixMaskQuickTest<double> dense_matrix_mask_quick_test_double("double");
