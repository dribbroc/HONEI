/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_row_sum_vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

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
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index() + 1);
                }

                std::tr1::shared_ptr<DenseVector<DataType_> > dv(MatrixRowSumVector<DataType_>::value(*dm));
                DataType_ s(size);
                for (typename Vector<DataType_>::ElementIterator v(dv->begin_elements()), v_end(dv->end_elements()) ;
                        v != v_end ; ++v)
                {
                    DataType_ last((v.index() + 1) * s), first(v.index() * s);
                    DataType_ ssum((last * (last + 1) / 2) - (first * (first + 1) / 2));

                    TEST_CHECK_EQUAL_WITHIN_EPS(*v, ssum, ssum * 100 * std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseMatrixRowSumVectorTest<float> dense_matrix_row_sum_vector_test_float("float");
DenseMatrixRowSumVectorTest<double> dense_matrix_row_sum_vector_test_double("double");
