/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/matrix_element_inverse.hh>
#include <unittest/unittest.hh>

#include <iostream>
#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseMatrixElementInverseTest :
    public BaseTest
{
    public:
        DenseMatrixElementInverseTest(const std::string & type) :
            BaseTest("matrix_element_inverse_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(size, size,
                            static_cast<DataType_>(0))),
                        dm2(new DenseMatrix<DataType_>(size, size, static_cast<DataType_>(0)));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1->begin_elements()), i_end(dm1->end_elements()),
                        j(dm2->begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }

                TEST_CHECK_EQUAL(MatrixElementInverse<DataType_>::value(*dm1), *dm2);
            }
        }
};

DenseMatrixElementInverseTest<float> dense_matrix_element_inverse_test_float("float");
DenseMatrixElementInverseTest<double> dense_matrix_element_inverse_test_double("double");

template <typename DataType_>
class DenseMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_matrix_element_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(3, 2,
                        static_cast<DataType_>(0))),
                    dm2(new DenseMatrix<DataType_>(3, 2, static_cast<DataType_>(0)));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1->begin_elements()), i_end(dm1->end_elements()),
                    j(dm2->begin_elements()) ; i != i_end ; ++i)
            {
                *i = i.index() + 1;
                *j = 1 / static_cast<DataType_>(i.index() + 1);
                ++j;
            }

            TEST_CHECK_EQUAL(MatrixElementInverse<DataType_>::value(*dm1), *dm2);
        }
};
DenseMatrixElementInverseQuickTest<float>  dense_matrix_element_inverse_quick_test_float("float");
DenseMatrixElementInverseQuickTest<double> dense_matrix_element_inverse_quick_test_double("double");
