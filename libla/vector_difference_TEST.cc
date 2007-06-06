/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_difference.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorDifferenceTest :
    public BaseTest
{
    public:
        DenseVectorDifferenceTest(const std::string & type) :
            BaseTest("dense_vector_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size));
                    
                for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                for (typename Vector<DataType_>::ElementIterator i(dv2->begin_elements()), i_end(dv2->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                DenseVector<DataType_> difference1(VectorDifference<DataType_>::value(*dv1, *dv2));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(difference1));
                TEST_CHECK_EQUAL(v1, 0);
            }

            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(5,
                    static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorDifference<DataType_>::value(*dv00, *dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorDifferenceTest<float> dense_vector_difference_test_float("float");
DenseVectorDifferenceTest<double> dense_vector_difference_test_double("double");
