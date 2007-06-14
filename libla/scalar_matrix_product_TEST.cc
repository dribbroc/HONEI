/* vim: set sw=4 sts=4 et foldmethod=syntax : */
///< \todo Fix compile errors

#include <libla/dense_matrix.hh>
#include <libla/scalar_matrix_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class ScalarMatrixProductTest :
    public BaseTest
{
    public:
        ScalarMatrixProductTest(const std::string & type) :
            BaseTest("scalar_matrix_product_test<" + type + ">")
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

ScalarMatrixProductTest<float> scalar_matrix_product_test_float("float");
ScalarMatrixProductTest<double> scalar_matrix_product_test_double("double");
