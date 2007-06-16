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
class DenseScalarVectorProductTest :
    public BaseTest
{
    public:
        DenseScalarVectorProductTest(const std::string & type) :
            BaseTest("dense_scalar_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(3)));

                DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
                TEST_CHECK_EQUAL(v1, 6 * size);
            }
        }
};

DenseScalarVectorProductTest<float> dense_scalar_vector_product_test_float("float");
DenseScalarVectorProductTest<double> dense_scalar_vector_product_test_double("double");

template <typename DataType_>
class DenseScalarVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseScalarVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_scalar_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                static_cast<DataType_>(3)));

            DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
            TEST_CHECK_EQUAL(v1, 6 * size);
        }
};
DenseScalarVectorProductQuickTest<float>  dense_scalar_vector_product_quick_test_float("float");
DenseScalarVectorProductQuickTest<double> dense_scalar_vector_product_quick_test_double("double");
