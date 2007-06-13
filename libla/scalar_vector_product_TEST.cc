/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <scalar_vector_product.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class ScalarVectorProductTest :
    public BaseTest
{
    public:
        ScalarVectorProductTest(const std::string & type) :
            BaseTest("scalar_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));

                DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
                TEST_CHECK_EQUAL(v1, 2 * size);
            }
        }
};

ScalarVectorProductTest<float> scalar_vector_product_test_float("float");
ScalarVectorProductTest<double> scalar_vector_product_test_double("double");

template <typename DataType_>
class ScalarVectorProductQuickTest :
    public QuickTest
{
    public:
        ScalarVectorProductQuickTest(const std::string & type) :
            QuickTest("scalar_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                static_cast<DataType_>(1)));

            DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
            TEST_CHECK_EQUAL(v1, 2 * size);
        }
};
ScalarVectorProductQuickTest<float>  scalar_vector_product_quick_test_float("float");
ScalarVectorProductQuickTest<double> scalar_vector_product_quick_test_double("double");
