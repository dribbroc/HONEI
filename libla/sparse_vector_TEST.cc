/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <string>
#include <tr1/memory>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class SparseVectorCreationTest :
    public BaseTest
{
    public:
        SparseVectorCreationTest(const std::string & type) :
            BaseTest("sparse_vector_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<SparseVector<DataType_> > sm1(new SparseVector<DataType_>(size, size));
                std::tr1::shared_ptr<SparseVector<DataType_> > sm2(new SparseVector<DataType_>(size, 1));
                TEST_CHECK(true);
            }
        }
};

SparseVectorCreationTest<float> sparse_vector_creation_test_float("float");
SparseVectorCreationTest<double> sparse_vector_creation_test_double("double");
