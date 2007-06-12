/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_norm.hh>
#include <libla/scalar_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseScalarProductQuickTest :
    public QuickTest
{
    public:
        DenseScalarProductQuickTest(const std::string & type) :
            QuickTest("dense_scalar_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            int quicksize = 2 << 7;
            std::tr1::shared_ptr<DenseVector<DataType_> > dv0q(new DenseVector<DataType_>(quicksize,
                static_cast<DataType_>(0)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1q(new DenseVector<DataType_>(quicksize,
                static_cast<DataType_>(1)));

            DataType_ p0q(ScalarProduct<DataType_>::value(*dv1q,*dv0q));
            DataType_ p1q(ScalarProduct<DataType_>::value(*dv1q,*dv1q));
            TEST_CHECK_EQUAL(p0q,0);
            TEST_CHECK_EQUAL(p1q,quicksize);
            
            std::tr1::shared_ptr<DenseVector<DataType_> > dv2q(new DenseVector<DataType_>(quicksize));
            for (typename Vector<DataType_>::ElementIterator iq(dv2q->begin_elements()), i_endq(dv2q->end_elements()) ;
                    iq != i_endq ; ++iq)
            {
                *iq = static_cast<DataType_>((iq.index() + 1) / 1.23456789);
            }

            DataType_ v2q(VectorNorm<DataType_, vnt_l_two, false>::value(*dv2q));
            DataType_ p2q(ScalarProduct<DataType_>::value(*dv2q,*dv2q));
            TEST_CHECK_EQUAL_WITHIN_EPS(v2q, p2q,sqrt(std::numeric_limits<DataType_>::epsilon()));
        }
};
DenseScalarProductQuickTest<float>  dense_scalar_product_quick_test_float("float");
DenseScalarProductQuickTest<double> dense_scalar_product_quick_test_double("double");


template <typename DataType_>
class DenseScalarProductTest :
    public BaseTest
{
    public:
        DenseScalarProductTest(const std::string & type) :
            BaseTest("dense_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv0(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(0)));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));

                DataType_ p0(ScalarProduct<DataType_>::value(*dv1,*dv0));
                DataType_ p1(ScalarProduct<DataType_>::value(*dv1,*dv1));
                TEST_CHECK_EQUAL(p0,0);
                TEST_CHECK_EQUAL(p1,size);
                
                std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>(size));
                for (typename Vector<DataType_>::ElementIterator i(dv2->begin_elements()), i_end(dv2->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v2(VectorNorm<DataType_, vnt_l_two, false>::value(*dv2));
                DataType_ p2(ScalarProduct<DataType_>::value(*dv2,*dv2));
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2,sqrt(std::numeric_limits<DataType_>::epsilon()));
            }
            
            std::tr1::shared_ptr<DenseVector<DataType_> > dv00(new DenseVector<DataType_>(1,
                    static_cast<DataType_>(1)));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv01(new DenseVector<DataType_>(2,
                    static_cast<DataType_>(1)));
                    
            TEST_CHECK_THROWS(ScalarProduct<DataType_>::value(*dv00,*dv01),VectorSizeDoesNotMatch);            
        }
};

DenseScalarProductTest<float> dense_scalar_product_test_float("float");
DenseScalarProductTest<double> dense_scalar_product_double("double");
