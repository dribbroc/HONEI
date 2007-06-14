/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libla/scalar_matrix_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class ScalarDenseMatrixProductTest :
    public BaseTest
{
    public:
        ScalarDenseMatrixProductTest(const std::string & type) :
            BaseTest("scalar_dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(1)), dm2(size, size + 1, DataType_(2));
                DenseMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(2), dm1));

                for (typename Matrix<DataType_>::ConstElementIterator i(dm1.begin_elements()),
                        i_end(dm1.end_elements()), j(dm2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

ScalarDenseMatrixProductTest<float> scalar_dense_matrix_product_test_float("float");
ScalarDenseMatrixProductTest<double> scalar_dense_matrix_product_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_dense_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm(size, size + 1, DataType_(1));
            DenseMatrix<DataType_> prod1(ScalarMatrixProduct<DataType_>::value(DataType_(2), dm));

            DataType_ vsum(0), ssum(2 * DataType_(size) * DataType_(size + 1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                vsum += *i;
            }
            TEST_CHECK_EQUAL_WITHIN_EPS(vsum, static_cast<DataType_>(2 * size * (size +1)), 
                std::numeric_limits<DataType_>::epsilon());
        }
};
//ScalarDenseMatrixProductQuickTest<float> scalar_dense_matrix_product_quick_test_float("float");
//ScalarDenseMatrixProductQuickTest<double> scalar_dense_matrix_product_quick_test_double("double");
