/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_absolute.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorAbsoluteValueTest :
    public BaseTest
{
    public:
        DenseVectorAbsoluteValueTest(const std::string & type) :
            BaseTest("dense_vector_absolute_value_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size)),
                        dv2(new DenseVector<DataType_>(size));
                for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                    (*dv2)[i.index()] = i.index() / 0.987654321;
                }

                TEST_CHECK_EQUAL(VectorAbsolute<DataType_>::value(*dv1), *dv2);
            }
        }
};

DenseVectorAbsoluteValueTest<float> dense_vector_absolute_value_test_float("float");
DenseVectorAbsoluteValueTest<double> dense_vector_absolute_value_test_double("double");

template <typename DataType_>
class DenseVectorAbsoluteQuickTest :
    public QuickTest
{
    public:
        DenseVectorAbsoluteQuickTest(const std::string & type) :
            QuickTest("dense_vector_absolute_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size)),
                    dv2(new DenseVector<DataType_>(size));
            for (typename Vector<DataType_>::ElementIterator i(dv1->begin_elements()), i_end(dv1->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                (*dv2)[i.index()] = i.index() / 0.987654321;
            }

            TEST_CHECK_EQUAL(VectorAbsolute<DataType_>::value(*dv1), *dv2);
        }
};
DenseVectorAbsoluteQuickTest<float>  dense_vector_absolute_quick_test_float("float");
DenseVectorAbsoluteQuickTest<double> dense_vector_absolute_quick_test_double("double");
