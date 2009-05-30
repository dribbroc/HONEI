/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <unittest/unittest.hh>

#include <string>
#include <tr1/memory>

using namespace honei;
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
            //for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                unsigned long size(5);
                unsigned long columns(2 * size), rows(size);
                SparseMatrix<DataType_> sm1(rows, columns, size / 10 + 1);
                SparseMatrix<DataType_> sm2(sm1.copy());

                TEST_CHECK_EQUAL(sm2.rows(), sm1.rows());
                TEST_CHECK_EQUAL(sm2.columns(), sm1.columns());

                for (typename SparseMatrix<DataType_>::ElementIterator i(sm2.begin_elements()),
                        i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    typename SparseMatrix<DataType_>::ConstElementIterator ci(i);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*ci, DataType_(0), std::numeric_limits<DataType_>::epsilon());

                    if (0 == (i.index() % 7))
                    {
                        *i = i.index();
                    }
                }

                for (typename SparseMatrix<DataType_>::ConstElementIterator i(sm1.begin_elements()),
                        i_end(sm1.end_elements()), j(sm2.begin_elements()) ; i != i_end ; ++i, ++j)
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
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);
                SparseMatrix<DataType_> sv1(rows, columns, columns);
                SparseMatrix<DataType_> sv2(1,size);
                TEST_CHECK(true);
            }
        }
};
SparseMatrixCreationTest<float> sparse_matrix_creation_test_float("float");
SparseMatrixCreationTest<double> sparse_matrix_creation_test_double("double");

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
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);

                SparseMatrix<DataType_> sm(rows, columns, size / 10 + 1);
                for (typename SparseMatrix<DataType_>::ElementIterator i(sm.begin_elements()), i_end(sm.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index());
                }

                TEST_CHECK_EQUAL(sm.columns(), columns);
                TEST_CHECK_EQUAL(sm.rows(), rows);

                SparseVector<DataType_> & row1 = sm[0];
                TEST_CHECK_EQUAL(row1.size(), columns);

                for (typename SparseVector<DataType_>::ConstElementIterator i(row1.begin_elements()), i_end(row1.end_elements()) ;
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
            SparseMatrix<DataType_> sm(rows, columns, 1);

            TEST_CHECK_EQUAL(sm.columns(), columns);
            TEST_CHECK_EQUAL(sm.rows(), rows);

            SparseVector<DataType_> & row1 = sm[0];

            TEST_CHECK_EQUAL(row1.size(), columns);

            for (typename SparseMatrix<DataType_>::NonZeroElementIterator i(sm.begin_non_zero_elements()),
                i_end(sm.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                //iterating over an empty matrix - should never reach this point
                TEST_CHECK(false);
            }


            unsigned long size (5);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1);
            for (typename SparseMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()); i != i_end ; ++i)
            {
                if (i.index() % 4 == 0) 
                {
                    *i = DataType_(5);

                }
            }
            for (typename SparseMatrix<DataType_>::NonZeroConstElementIterator i(sm1.begin_non_zero_elements()),
                i_end(sm1.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL(i.index() % 4, 0ul);
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(5), std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL(i_end.index(), size * (size + 1));
            }


            SparseMatrix<DataType_> sm3(size, size + 4);
            sm3(0, 1) = 1;
            sm3(1, 2) = 2;
            sm3(2, 1) = 3;
            sm3(2, 3) = 4;
            sm3(4, 4) = 5;
            DenseMatrix<DataType_> dm0(sm3);
            SparseMatrix<DataType_> sm4(dm0);
            TEST_CHECK_EQUAL(sm4, sm3);
        }
};
SparseMatrixQuickTest<float>  sparse_matrix_quick_test_float("float");
SparseMatrixQuickTest<double> sparse_matrix_quick_test_double("double");
