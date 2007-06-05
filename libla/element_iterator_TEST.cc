/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/element_iterator.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class BandedMatrixElementIterationTest :
    public BaseTest
{
    public:
        BandedMatrixElementIterationTest(const std::string & type) :
            BaseTest("banded_matrix_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<BandedMatrix<DataType_> > bm(new BandedMatrix<DataType_>(size));

                typename Matrix<DataType_>::ConstElementIterator ce(bm->begin_elements()), ce_end(bm->end_elements());
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 0, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

BandedMatrixElementIterationTest<float> banded_matrix_element_iteration_test_float("float");
BandedMatrixElementIterationTest<double> banded_matrix_element_iteration_test_double("double");

template <typename DataType_>
class DenseMatrixElementIterationTest :
    public BaseTest
{
    public:
        DenseMatrixElementIterationTest(const std::string & type) :
            BaseTest("dense_matrix_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(size, size,
                            static_cast<DataType_>(10)));

                typename MutableMatrix<DataType_>::ElementIterator e(dm->begin_elements()), e_end(dm->end_elements()) ;
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(e.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*e, 10, std::numeric_limits<DataType_>::epsilon());
                    *e = 333;
                    ++e;
                }

                for (typename MutableMatrix<DataType_>::ElementIterator me(dm->begin_elements()), me_end(dm->end_elements()) ;
                        me != me_end ; ++me)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*me, 333, std::numeric_limits<DataType_>::epsilon());
                }

                typename Matrix<DataType_>::ConstElementIterator ce(dm->begin_elements()), ce_end(dm->end_elements()) ;
                for (unsigned long i(0) ; i < (size * size) ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 333, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

DenseMatrixElementIterationTest<float> dense_matrix_element_iteration_test_float("float");
DenseMatrixElementIterationTest<double> dense_matrix_element_iteration_test_double("double");

template <typename DataType_>
class DenseVectorElementIterationTest :
    public BaseTest
{
    public:
        DenseVectorElementIterationTest(const std::string & type) :
            BaseTest("dense_vector_element_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size, static_cast<DataType_>(10)));

                typename Vector<DataType_>::ElementIterator e(dv->begin_elements()), e_end(dv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(e.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*e, 10, std::numeric_limits<DataType_>::epsilon());
                    *e = 222;
                    ++e;
                }

                typename Vector<DataType_>::ConstElementIterator ce(dv->begin_elements()), ce_end(dv->end_elements()) ;
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL(ce.index(), i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 222, std::numeric_limits<DataType_>::epsilon());
                    ++ce;
                }
            }
        }
};

DenseVectorElementIterationTest<float> dense_vector_element_iteration_test_float("float");
DenseVectorElementIterationTest<double> dense_vector_element_iteration_test_double("double");
