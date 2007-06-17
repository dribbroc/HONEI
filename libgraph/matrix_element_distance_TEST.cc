/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libgraph/matrix_element_distance.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;
template <typename DataType_>
class MatrixElementDistanceQuickTest :
    public QuickTest
{
    public:
        MatrixElementDistanceQuickTest(const std::string & type) :
            QuickTest("matrix_element_distance_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long columns(3), rows(2);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>
                (3, 2));

            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()),
                    i_end(dm->end_elements()) ; i != i_end ; ++i)
            {
                *i = i.index();
            }
            DenseMatrix<DataType_> distance = MatrixElementDistance::value(*dm);
            
            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);            
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][1], sqrt(8), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][0], sqrt(8), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][2], sqrt(32), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][2], sqrt(8), std::numeric_limits<DataType_>::epsilon());  
            
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm2(new DenseMatrix<DataType_>
                (3, 3));            
            TEST_CHECK_THROWS(MatrixElementDIstance::value(*dm2), MatrixRowsDoNotMatch);
                                                                                              
        }
};
MatrixElementDistanceQuickTest<float>  matrix_element_distance_quick_test_float("float");
MatrixElementDistanceQuickTest<double> matrix_element_distance_test_quick_double("double");
