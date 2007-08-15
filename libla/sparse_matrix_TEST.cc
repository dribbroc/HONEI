/* vim: set sw=4 sts=4 et folsmethod=syntax : */

#include <libla/sparse_matrix.hh>
#include <unittest/unittest.hh>

#include <string>
#include <tr1/memory>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class SparseMatrixCopyTest :
    public BaseTest
{
    public:
        SparseMatrixCopyTest(const std::string & type) :
            BaseTest("sparse_matrix_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                unsigned long columns(2 * size), rows(size);
                SparseMatrix<DataType_> sm1(columns, rows, size / 10 + 1);
                std::tr1::shared_ptr<SparseMatrix<DataType_> > sm2(sm1.copy());

                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2->begin_elements()),
                        i_end(sm2->end_elements()) ; i != i_end ; ++i)
                {
                    typename Matrix<DataType_>::ConstElementIterator ci(i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ci, DataType_(0), std::numeric_limits<DataType_>::epsilon());

                    if (0 == (i.index() % 7))
                        *i = i.index();
                }

                for (typename Matrix<DataType_>::ConstElementIterator i(sm1.begin_elements()),
                        i_end(sm1.end_elements()), j(sm2->begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());

                    DataType_ s(0 == (j.index() % 7) ? j.index() : 0);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*j, s, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
SparseMatrixCopyTest<float> sparse_matrix_copy_test_float("float");
SparseMatrixCopyTest<double> sparse_matrix_copy_test_double("double");

template <typename DataType_>
class SparseMatrixCreationTest :
    public BaseTest
{
    public:
        SparseMatrixCreationTest(const std::string & type) :
            BaseTest("sparse_matrix_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);
                SparseMatrix<DataType_> sv1(columns, rows, columns * rows);
                SparseMatrix<DataType_> sv2(size, 1);
                TEST_CHECK(true);
            }
        }
};
SparseMatrixCreationTest<float> sparse_matrix_creation_test_float("float");
SparseMatrixCreationTest<double> sparse_matrix_creation_test_double("double");
SparseMatrixCreationTest<int> sparse_matrix_creation_test_int("int");
SparseMatrixCreationTest<bool> sparse_matrix_creation_test_bool("bool");

template <typename DataType_>
class SparseMatrixLayoutTest :
    public BaseTest
{
    public:
        SparseMatrixLayoutTest(const std::string & type) :
            BaseTest("sparse_matrix_layout_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);

                SparseMatrix<DataType_> sm(columns, rows, size / 10 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm.begin_elements()), i_end(sm.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index());
                }

                TEST_CHECK_EQUAL(sm.columns(), columns);
                TEST_CHECK_EQUAL(sm.rows(), rows);

                Vector<DataType_> & row1 = sm[0];
                TEST_CHECK_EQUAL(row1.size(), columns);

                for (typename Vector<DataType_>::ConstElementIterator i(row1.begin_elements()), i_end(row1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, i.index(), std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

SparseMatrixLayoutTest<float> sparse_matrix_layout_test_float("float");
SparseMatrixLayoutTest<double> sparse_matrix_layout_test_double("double");

template <typename DataType_>
class SparseMatrixQuickTest :
    public QuickTest
{
    public:
        SparseMatrixQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long columns(3), rows(2);
            SparseMatrix<DataType_> sm(columns, rows, 1);

            TEST_CHECK_EQUAL(sm.columns(), columns);
            TEST_CHECK_EQUAL(sm.rows(), rows);

            Vector<DataType_> & row1 = sm[0];

            TEST_CHECK_EQUAL(row1.size(), columns);

            for (typename MutableMatrix<DataType_>::ElementIterator i(sm.begin_non_zero_elements()),
                i_end(sm.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                //iterating over an empty matrix - should never reach this point
                TEST_CHECK(false);
            } 
            
            sm[0][0] = DataType_(1);
            sm[0][1] = DataType_(2);            
            sm[1][1] = DataType_(3);
            typename MutableMatrix<DataType_>::ElementIterator i(sm.begin_non_zero_elements());
            typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_non_zero_elements());            
            TEST_CHECK_EQUAL(*i, DataType_(1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm.begin_non_zero_elements()),
                i_end(sm.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL(*i, i.index()+1);
            }            
            
        }
};
SparseMatrixQuickTest<float>  sparse_matrix_quick_test_float("float");
SparseMatrixQuickTest<double> sparse_matrix_quick_test_double("double");
