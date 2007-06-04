/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseMatrixCreationTest :
    public BaseTest
{
    public:
        DenseMatrixCreationTest(const std::string & type) :
            BaseTest("dense_matrix_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size, 
                            static_cast<DataType_>(0)));
                TEST_CHECK(true);
            }
        }
};

DenseMatrixCreationTest<float> dense_matrix_creation_test_float("float");
DenseMatrixCreationTest<double> dense_matrix_creation_test_double("double");
