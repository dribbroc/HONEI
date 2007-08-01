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
#include <libla/matrix_row_sum_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class BandedMatrixRowSumVectorTest :
    public BaseTest
{
    public:
        BandedMatrixRowSumVectorTest(const std::string & type) :
            BaseTest("banded_matrix_row_sum_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                    
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseVector<DataType_> sum(MatrixRowSumVector<DataType_>::value(bm1));
                DenseVector<DataType_> dv2(size, static_cast<DataType_>(2));
                
                TEST_CHECK_EQUAL(sum, dv2);
            }
        }
};
BandedMatrixRowSumVectorTest<float> banded_matrix_row_sum_vector_test_float("float");
BandedMatrixRowSumVectorTest<double> banded_matrix_row_sum_vector_test_double("double");

template <typename DataType_>
class BandedMatrixRowSumVectorQuickTest :
    public QuickTest
{
    public:
        BandedMatrixRowSumVectorQuickTest(const std::string & type) :
            QuickTest("banded_matrix_row_sum_vector_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                    
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseVector<DataType_> sum(MatrixRowSumVector<DataType_>::value(bm1));
            DenseVector<DataType_> dv2(size, static_cast<DataType_>(2));
            
            TEST_CHECK_EQUAL(sum, dv2);
        }
};
BandedMatrixRowSumVectorQuickTest<float> banded_matrix_row_sum_vector_quick_test_float("float");
BandedMatrixRowSumVectorQuickTest<double> banded_matrix_row_sum_vector_quick_test_double("double");

template <typename DataType_>
class DenseMatrixRowSumVectorTest :
    public BaseTest
{
    public:
        DenseMatrixRowSumVectorTest(const std::string & type) :
            BaseTest("dense_matrix_row_sum_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index() + 1);
                }

                DenseVector<DataType_> dv(MatrixRowSumVector<DataType_>::value(*dm));
                DataType_ s(size);
                for (typename Vector<DataType_>::ElementIterator v(dv.begin_elements()), v_end(dv.end_elements()) ;
                        v != v_end ; ++v)
                {
                    DataType_ last((v.index() + 1) * s), first(v.index() * s);
                    DataType_ ssum((last * (last + 1) / 2) - (first * (first + 1) / 2));

                    TEST_CHECK_EQUAL_WITHIN_EPS(*v, ssum, ssum * 200 * std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseMatrixRowSumVectorTest<float> dense_matrix_row_sum_vector_test_float("float");
DenseMatrixRowSumVectorTest<double> dense_matrix_row_sum_vector_test_double("double");

template <typename DataType_>
class DenseMatrixRowSumVectorQuickTest :
    public QuickTest
{
    public:
        DenseMatrixRowSumVectorQuickTest(const std::string & type) :
            QuickTest("dense_matrix_row_sum_vector_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(3);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size+2));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>(i.index() + 1);
            }

            DenseVector<DataType_> dv(MatrixRowSumVector<DataType_>::value(*dm));
            DataType_ s(size);
            for (typename Vector<DataType_>::ElementIterator v(dv.begin_elements()), v_end(dv.end_elements()) ;
                    v != v_end ; ++v)
            {
                DataType_ last((v.index() + 1) * s), first(v.index() * s);
                DataType_ ssum((last * (last + 1) / 2) - (first * (first + 1) / 2));

                TEST_CHECK_EQUAL_WITHIN_EPS(*v, ssum, ssum * 100 * std::numeric_limits<DataType_>::epsilon());
            }
        }
};
DenseMatrixRowSumVectorQuickTest<float>  dense_matrix_row_sum_vector_quick_test_float("float");
DenseMatrixRowSumVectorQuickTest<double> dense_matrix_row_sum_vector_quick_test_double("double");

template <typename DataType_>
class SparseMatrixRowSumVectorTest :
    public BaseTest
{
    public:
        SparseMatrixRowSumVectorTest(const std::string & type) :
            BaseTest("sparse_matrix_row_sum_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 10) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % (size) == 0) 
                    {
                        *i = DataType_(2);
                    }
                }   
                DenseVector<DataType_> sum(MatrixRowSumVector<DataType_>::value(sm1));
                DenseVector<DataType_> dv1(size + 1, static_cast<DataType_>((size / size / 2 + 1) * 2));
                
                TEST_CHECK_EQUAL(sum, dv1);
            }
        }
};
SparseMatrixRowSumVectorTest<float> sparse_matrix_row_sum_vector_test_float("float");
SparseMatrixRowSumVectorTest<double> sparse_matrix_row_sum_vector_test_double("double");

template <typename DataType_>
class SparseMatrixRowSumVectorQuickTest :
    public QuickTest
{
    public:
        SparseMatrixRowSumVectorQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_row_sum_vector_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % ((size + 1) / 2) == 0) 
                {
                    *i = DataType_(2);
                }
            }   
            DenseVector<DataType_> sum(MatrixRowSumVector<DataType_>::value(sm1));
            DenseVector<DataType_> dv1(size + 1, static_cast<DataType_>((size / size / 2 + 1) * 4));
            
            TEST_CHECK_EQUAL(sum, dv1);
        }
};
SparseMatrixRowSumVectorQuickTest<float> sparse_matrix_row_sum_vector_quick_test_float("float");
SparseMatrixRowSumVectorQuickTest<double> sparse_matrix_row_sum_vector_quick_test_double("double");
