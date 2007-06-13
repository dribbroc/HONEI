/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
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
class MatrixMastTest :
    public BaseTest
{
    public:
        MatrixMastTest(const std::string & type) :
            BaseTest("matrix_mask_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                std::tr1::shared_ptr<DenseMatrix<DataType_> > mask(new DenseMatrix<bool>(
                    size, size+1, false));
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(
                    size, size+1, static_cast<DataType_>(1)));
                
                for (typename MutableMatrix<DataType_>::ElementIterator i(mask->begin_elements()), 
                    i_end(mask->end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index()%2 == 0) 
                    {
                        *i=true;
                        ++count;
                    }
                    
                }

                DenseMatrix<DataType_> result(MatrixMask<DataType_>::value(*mask , *dm));
                DataType_ sum = 0;
                for (typename MutableMatrix<DataType_>::ElementIterator j(result.begin_elements()), 
                    j_end(result.end_elements()) ; j != j_end ; ++j)
                {
                    sum += *j;
                }
                
                TEST_CHECK_EQUAL(sum, count);
            }
            
            std::tr1::shared_ptr<Matrix<DataType_> > mask1(new DenseMatrix<bool>(
                2, 3, false));
            std::tr1::shared_ptr<Matrix<DataType_> > mask2(new DenseMatrix<bool>(
                3, 4, false));                
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(
                2, 4, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask1), MatrixRowsDoesNotMatch);
            TEST_CHECK_THROWS(MatrixMask<DataType_>::value(*dm1, *mask2), MatrixColumnsDoesNotMatch);
        }
};

MatrixMastTest<float> matrix_mask_test_float("float");
MatrixMastTest<double> matrix_mask_test_double("double"); */
