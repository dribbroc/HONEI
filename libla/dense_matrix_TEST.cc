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
class DenseMatrixCopyTest :
    public BaseTest
{
    public:
        DenseMatrixCopyTest(const std::string & type) :
            BaseTest("dense_vector_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 4) ; size <<= 1)
            {
                DenseMatrix<DataType_> dv1(size, size, static_cast<DataType_>(0)),
                    dv2(size, size, static_cast<DataType_>(1));
                std::tr1::shared_ptr<DenseMatrix<DataType_> > c(dv1.copy());

                for (typename MutableMatrix<DataType_>::ElementIterator i(c->begin_elements()),
                        i_end(c->end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename Matrix<DataType_>::ConstElementIterator i(dv1.begin_elements()),
                        i_end(dv1.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }

                for (typename Matrix<DataType_>::ConstElementIterator i(dv2.begin_elements()),
                        i_end(dv2.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 1, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseMatrixCopyTest<float> dense_matrix_copy_test_float("float");
DenseMatrixCopyTest<double> dense_matrix_copy_test_double("double");

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

                TEST_CHECK_EQUAL(*dm0, *dm1);
            }
        }
};
DenseMatrixEqualityTest<bool> dense_matrix_equality_test_float("float");
DenseMatrixEqualityTest<int> dense_matrix_equality_test_double("double");

template <typename DataType_>
class DenseMatrixLayoutTest :
    public BaseTest
{
    public:
        DenseMatrixLayoutTest(const std::string & type) :
            BaseTest("dense_matrix_layout_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);

                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(columns, rows));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()), i_end(dm->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index());
                }


                TEST_CHECK_EQUAL(dm->columns(), columns);
                TEST_CHECK_EQUAL(dm->rows(), rows);

                Vector<DataType_> & row1 = (*dm)[0];
                Vector<DataType_> & col1 = dm->column(0);

                TEST_CHECK_EQUAL(row1.size(), columns);
                TEST_CHECK_EQUAL(col1.size(), rows);

                for (typename Vector<DataType_>::ConstElementIterator i(row1.begin_elements()), i_end(row1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, i.index(), std::numeric_limits<DataType_>::epsilon());
                }

                for (typename Vector<DataType_>::ConstElementIterator i(col1.begin_elements()), i_end(col1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, i.index() * columns, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseMatrixLayoutTest<float> dense_matrix_layout_test_float("float");
DenseMatrixLayoutTest<double> dense_matrix_layout_test_double("double");

template <typename DataType_>
class DenseMatrixQuickTest :
    public QuickTest
{
    public:
        DenseMatrixQuickTest(const std::string & type) :
            QuickTest("dense_matrix_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long columns(3), rows(2);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>
                (columns, rows,static_cast<DataType_>(1)));

            TEST_CHECK_EQUAL(*dm, *dm);
            TEST_CHECK_EQUAL(dm->columns(), columns);
            TEST_CHECK_EQUAL(dm->rows(), rows);

            Vector<DataType_> & row1 = (*dm)[0];
            Vector<DataType_> & col1 = dm->column(0);

            TEST_CHECK_EQUAL(row1.size(), columns);
            TEST_CHECK_EQUAL(col1.size(), rows);

        }
};
DenseMatrixQuickTest<float>  dense_matrix_quick_test_float("float");
DenseMatrixQuickTest<double> dense_matrix_quick_test_double("double");
