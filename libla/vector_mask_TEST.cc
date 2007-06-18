/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <libla/vector_mask.hh>
#include <libla/vector_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorMaskTest :
    public BaseTest
{
    public:
        DenseVectorMaskTest(const std::string & type) :
            BaseTest("dense_vector_mask_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ count = 0;
                std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                    size, false));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(
                    size, static_cast<DataType_>(1)));
                
                for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                    i_end(mask->end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index()%2 == 0) 
                    {
                        *i=true;
                        ++count;
                    }
                    
                }
                DenseVector<DataType_> result(VectorMask<DataType_>::value(*dv, *mask));
                DataType_ sum = 0;
                for (typename Vector<DataType_>::ElementIterator j(result.begin_elements()), 
                    j_end(result.end_elements()) ; j != j_end ; ++j)
                {
                    sum += *j;
                }
                
                TEST_CHECK_EQUAL(sum, count);
            }

            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMask<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskTest<float> dense_vector_mask_test_float("float");
DenseVectorMaskTest<double> dense_vector_mask_test_double("double"); 

template <typename DataType_>
class DenseVectorMaskQuickTest :
    public QuickTest
{
    public:
        DenseVectorMaskQuickTest(const std::string & type) :
            QuickTest("dense_vector_mask_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DataType_ count = 0;
            std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                size, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(
                size, static_cast<DataType_>(1)));
            
            for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                i_end(mask->end_elements()) ; i != i_end ; ++i)
            {
                if (i.index()%2 == 0) 
                {
                    *i=true;
                    ++count;
                }
                
            }
            DenseVector<DataType_> result(VectorMask<DataType_>::value(*dv, *mask));
            DataType_ sum = 0;
            for (typename Vector<DataType_>::ElementIterator j(result.begin_elements()), 
                j_end(result.end_elements()) ; j != j_end ; ++j)
            {
                sum += *j;
            }
            
            TEST_CHECK_EQUAL(sum, count);

            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMask<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);

        }
};
DenseVectorMaskQuickTest<float>  dense_vector_mask_quick_test_float("float");
DenseVectorMaskQuickTest<double> dense_vector_mask_quick_test_double("double");
