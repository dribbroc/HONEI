/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_row_sum_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;

template <typename DataType_>
class DenseMatrixRowSumVectorTest :
    public BaseTest
{
    public:
        DenseMatrixRowSumVectorTest(const std::string & type) :
            BaseTest("dense_matrix_row_sum_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size, 
                            static_cast<DataType_>(1)));
                
                DenseVector<DataType_> dv(MatrixRowSumVector<DataType_>::value(*dm));
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), i_end(dv->end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL(i.size(),1);
                }                
                            
            }
        }
};

DenseMatrixRowSumVectorTest<float> dense_matrix_row_sum_vector_test_float("float");
DenseMatrixRowSumVectorTest<double> dense_matrix_row_sum_vector_test_double("double");
