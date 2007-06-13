/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/banded_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

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
                std::tr1::shared_ptr<BandedMatrix<DataType_> > dm1(new BandedMatrix<DataType_>(size));
                TEST_CHECK(true);
                //std::tr1::shared_ptr<DenseVector<DataType_> >dv1(new DenseVector<DataType_>
                //    (size, static_cast<DataType_>(1)));
                //std::tr1::shared_ptr<BandedMatrix<DataType_> > dm2(new BandedMatrix<DataType_>(size,*dv1));
                //TEST_CHECK(true);
                
            }
            
            //std::tr1::shared_ptr<DenseVector<DataType_> > dv2(new DenseVector<DataType_>
            //    (5, static_cast<DataType_>(1)));                
            //TEST_CHECK_THROWS (std::tr1::shared_ptr<BandedMatrix<DataType_> > dm3(new BandedMatrix<DataType_>(6,*dv2)),
            //    VectorSizeDoesNotMatch);
        }
};

BandedMatrixCreationTest<float> banded_matrix_creation_test_float("float");
BandedMatrixCreationTest<double> banded_matrix_creation_test_double("double");

template <typename DataType_>
class BandedMatrixQuickTest :
    public QuickTest
{
    public:
        BandedMatrixQuickTest(const std::string & type) :
            QuickTest("banded_matrix_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<BandedMatrix<DataType_> > dv(new BandedMatrix<DataType_>(4711));
            TEST_CHECK(true);
        }
};
BandedMatrixQuickTest<float>  banded_matrix_quick_test_float("float");
BandedMatrixQuickTest<double> banded_matrix_quick_test_double("double");
