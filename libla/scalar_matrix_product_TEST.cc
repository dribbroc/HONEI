/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libla/scalar_matrix_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseScalarMatrixProductTest :
    public BaseTest
{
    public:
        DenseScalarMatrixProductTest(const std::string & type) :
            BaseTest("dense_scalar_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(
                    size, size+1, static_cast<DataType_>(1)));

                DenseMatrix<DataType_> prod1(ScalarMatrixProduct<DataType_>::value(static_cast<DataType_>(2), *dm));
                DataType_ sum = 0;
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                        i != i_end ; ++i)
                {
                    sum += *i;
                }
                
                TEST_CHECK_EQUAL(sum, static_cast<DataType_>(2 * size * (size +1)));
            }
        }
};

DenseScalarMatrixProductTest<float> dense_scalar_matrix_product_test_float("float");
DenseScalarMatrixProductTest<double> dense_scalar_matrix_product_test_double("double");

template <typename DataType_>
class DenseScalarMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseScalarMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_scalar_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(
                size, size+1, static_cast<DataType_>(1)));

            DenseMatrix<DataType_> prod1(ScalarMatrixProduct<DataType_>::value(static_cast<DataType_>(2), *dm));
            DataType_ sum = 0;
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                    i != i_end ; ++i)
            {
                sum += *i;
            }
            
            TEST_CHECK_EQUAL(sum, static_cast<DataType_>(2 * size * (size +1)));
        }
};
DenseScalarMatrixProductQuickTest<float>  dense_scalar_matrix_product_quick_test_float("float");
DenseScalarMatrixProductQuickTest<double> dense_scalar_matrix_product_quick_test_double("double");
