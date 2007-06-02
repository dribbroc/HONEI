/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/banded_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;

template <typename DataType_>
class BandedMatrixCreationTest :
    public BaseTest
{
    public:
        BandedMatrixCreationTest(const std::string & type) :
            BaseTest("banded_matrix_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<BandedMatrix<DataType_> > dm(new BandedMatrix<DataType_>(size));
                TEST_CHECK(true);
            }
        }
};

BandedMatrixCreationTest<float> banded_matrix_creation_test_float("float");
BandedMatrixCreationTest<double> banded_matrix_creation_test_double("double");
