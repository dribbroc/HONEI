/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;

template <typename DataType_>
class VectorNormValueTest :
    public BaseTest
{
    public:
        VectorNormValueTest(const std::string & type) :
            BaseTest("vector_norm_value_test<" + type + ">")
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

                DataType_ vmax(VectorNorm<DataType_, vnt_max>::value(*dv));
                DataType_ smax(size / 1.23456789);
                TEST_CHECK_EQUAL_WITHIN_EPS(vmax, smax, std::numeric_limits<DataType_>::epsilon());

                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(*dv));
                DataType_ s1(size * (size / 1.23456789 + 1) / 2);
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, std::numeric_limits<DataType_>::epsilon());

                DataType_ v2(VectorNorm<DataType_, vnt_l_two, false>::value(*dv));
                DataType_ s2(size * (size / 1.23456789 + 1) * (2 * size / 1.23456789 + 1) / 6);
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, s2, sqrt(std::numeric_limits<DataType_>::epsilon()));
            }
        }
};

VectorNormValueTest<float> vector_norm_value_test_float("float");
VectorNormValueTest<double> vector_norm_value_test_double("double");
