/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_masked_min_index.hh>
#include <libla/vector_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorMaskedMinIndexTest :
    public BaseTest
{
    public:
        DenseVectorMaskedMinIndexTest(const std::string & type) :
            BaseTest("dense_vector_masked_min_index_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                    size, false));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
                
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), 
                    i_end(dv->end_elements()) ; i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                
                (*dv)[3] = static_cast<DataType_>(0.000009);
                                    
                for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                    i_end(mask->end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index()%2 == 0) 
                    {
                        *i=true;
                        ++count;
                    }
                    
                }
                DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
                TEST_CHECK_EQUAL(result, (*dv)[3]);
                (*dv)[3] = static_cast<DataType_>(0);
                DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
                TEST_CHECK_EQUAL(result, (*dv)[3]);                
                (*dv)[3] = static_cast<DataType_>(-5.2345);
                DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
                TEST_CHECK_EQUAL(result, (*dv)[3]);                
            }

            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMinIndexTest<float> dense_vector_masked_min_index_test_float("float");
DenseVectorMaskedMinIndexTest<double> dense_vector_masked_min_index_test_double("double"); 

template <typename DataType_>
class DenseVectorMaskedMinIndexQuickTest :
    public QuickTest
{
    public:
        DenseVectorMaskedMinIndexQuickTest(const std::string & type) :
            QuickTest("dense_vector_masked_min_index_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DataType_ count = 0;
            std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                size, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
            
            for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), 
                i_end(dv->end_elements()) ; i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }
            
            (*dv)[3] = static_cast<DataType_>(0.000009);
                                
            for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                i_end(mask->end_elements()) ; i != i_end ; ++i)
            {
                if (i.index()%2 == 0) 
                {
                    *i=true;
                    ++count;
                }
                
            }
            DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
            TEST_CHECK_EQUAL(result, (*dv)[3]);
            (*dv)[3] = static_cast<DataType_>(0);
            DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
            TEST_CHECK_EQUAL(result, (*dv)[3]);                
            (*dv)[3] = static_cast<DataType_>(-5.2345);
            DataType_ result(VectorMaskedMinIndex<DataType_>::value(*dv, *mask));
            TEST_CHECK_EQUAL(result, (*dv)[3]);            


            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMaskedMinIndex<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMinIndexQuickTest<float>  dense_vector_masked_min_index_quick_test_float("float");
DenseVectorMaskedMinIndexQuickTest<double> dense_vector_masked_min_index_quick_test_double("double");
