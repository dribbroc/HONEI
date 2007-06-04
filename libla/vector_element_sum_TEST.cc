/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_element_sum.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <iostream>
#include <vector>
#include <tr1/memory>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorElementSumTest :
    public BaseTest
{
    public:
        DenseVectorElementSumTest(const std::string & type) :
            BaseTest("dense_vector_element_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), i_end(dv->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v1(VectorElementSum<DataType_>::value(*dv));
                DataType_ s1(size * (size + 1) / 2 / 1.23456789);
                // Behavious similar to size^2 * eps
                DataType_ eps1(s1 * 10 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

DenseVectorElementSumTest<float> dense_vector_element_sum_test_float("float");
DenseVectorElementSumTest<double> dense_vector_element_sum_test_double("double");
