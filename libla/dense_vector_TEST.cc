/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorCreationTest :
    public BaseTest
{
    public:
        DenseVectorCreationTest(const std::string & type) :
            BaseTest("dense_vector_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size,
                            static_cast<DataType_>(0)));
                TEST_CHECK(true);
            }
        }
};

DenseVectorCreationTest<float> dense_vector_creation_test_float("float");
DenseVectorCreationTest<double> dense_vector_creation_test_double("double");
DenseVectorCreationTest<bool> dense_vector_creation_test_bool("bool");
DenseVectorCreationTest<int> dense_vector_creation_test_int("int");

template <typename DataType_>
class DenseVectorEqualityTest :
    public BaseTest
{
    public:
        DenseVectorEqualityTest(const std::string & type) :
            BaseTest("dense_vector_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv0(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(1)));

                TEST_CHECK(dv0==dv1);
            }
        }
};
DenseVectorEqualityTest<bool> dense_vector_equality_test_bool("bool");
DenseVectorEqualityTest<int> dense_vector_equality_test_int("int");

template <typename DataType_>
class DenseVectorFunctionsTest :
    public BaseTest
{
    public:
        DenseVectorFunctionsTest(const std::string & type) :
            BaseTest("dense_vector_functions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), i_end(dv->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                TEST_CHECK_EQUAL(dv->size(),size);
                
                for (int i=0 ; i<size ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS( (*dv)[i],(i+1)/1.23456789, 
                        std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseVectorFunctionsTest<float> dense_vector_functions_test_float("float");
DenseVectorFunctionsTest<double> dense_vector_functions_test_double("double");

