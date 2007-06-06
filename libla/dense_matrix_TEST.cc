/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libla/vector.hh>
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

template <typename DataType_>
class DenseMatrixEqualityTest :
    public BaseTest
{
    public:
        DenseMatrixEqualityTest(const std::string & type) :
            BaseTest("dense_matrix_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm0(new DenseMatrix<DataType_>(size,
                    size,static_cast<DataType_>(1)));
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(size,
                    size,static_cast<DataType_>(1)));

                TEST_CHECK(dm0==dm1);
            }
        }
};
DenseMatrixEqualityTest<bool> dense_matrix_equality_test_bool("bool");
DenseMatrixEqualityTest<int> dense_matrix_equality_test_int("int");

template <typename DataType_>
class DenseMatrixFunctionsTest :
    public BaseTest
{
    public:
        DenseMatrixFunctionsTest(const std::string & type) :
            BaseTest("dense_matrix_functions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size,size+1));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), 
                    i_end(dm->end_elements()) ; i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index());
                }
                TEST_CHECK_EQUAL(dm->columns(),size);
                TEST_CHECK_EQUAL(dm->rows(),size+1);

                Vector<DataType_> * row_vector = &((*dm)[0]);
                Vector<DataType_> * column_vector = &(dm->column(0));                
                TEST_CHECK_EQUAL(row_vector->size(),size+1);
                TEST_CHECK_EQUAL(column_vector->size(),size);
                
                ///< \todo Fix compile errors in vector.hh and expand check over all elements of dm
                //TEST_CHECK_EQUAL_WITHIN_EPS(row_vector[1],1,std::numeric_limits<DataType_>::epsilon());

            }
        }
};

DenseMatrixFunctionsTest<float> dense_matrix_functions_test_float("float");
DenseMatrixFunctionsTest<double> dense_matrix_functions_test_double("double");
