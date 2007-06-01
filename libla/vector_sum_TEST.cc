/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_sum.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;

template <typename DataType_>
class DenseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sum_test<" + type + ">")
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

                DataType_ s1(VectorSum<DataType_>::value(*dv));
                DataType_ sum(size * (size / 1.23456789 + 1) / 2);
                TEST_CHECK_EQUAL_WITHIN_EPS(s1, sum, std::numeric_limits<DataType_>::epsilon());
            }
        }
};

DenseVectorSumTest<float> dense_vector_sum_test_float("float");
DenseVectorSumTest<double> dense_vector_sum_test_double("double");
