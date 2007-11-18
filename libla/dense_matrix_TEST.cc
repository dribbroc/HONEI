/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/vector.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
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
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm(size, size, DataType_(0));
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
                DenseMatrix<DataType_> dm1(size+1, size, DataType_(0)),
                    dm2(size+1, size, static_cast<DataType_>(1));
                std::tr1::shared_ptr<DenseMatrix<DataType_> > c(dm1.copy());

                for (typename MutableMatrix<DataType_>::ElementIterator i(c->begin_elements()),
                        i_end(c->end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename Matrix<DataType_>::ConstElementIterator i(dm1.begin_elements()),
                        i_end(dm1.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }

                for (typename Matrix<DataType_>::ConstElementIterator i(dm2.begin_elements()),
                        i_end(dm2.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 1, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseMatrixCopyTest<float> dense_matrix_copy_test_float("float");
DenseMatrixCopyTest<double> dense_matrix_copy_test_double("double");

template <typename DataType_>
class DenseMatrixDensifyQuickTest :
    public QuickTest
{
    public:
        DenseMatrixDensifyQuickTest(const std::string & type) :
            QuickTest("dense_matrix_densify_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
                unsigned long size(47);
                DenseMatrix<DataType_> dm0(size + 2, size, DataType_(0));
                SparseMatrix<DataType_> sm0(size + 2, size, size / 8 + 1);

                for (typename MutableMatrix<DataType_>::ElementIterator i(dm0.begin_elements()), j(sm0.begin_elements()),
                    i_end(dm0.end_elements()) ; i != i_end ; ++i , ++j)
                {
                    if (i.index() % 7 == 0)
                    {
                        *i = i.index();
                        *j = i.index();
                    }
                }
                DenseMatrix<DataType_> dm1(sm0);

                TEST_CHECK_EQUAL(dm1, dm0);

            }
};
DenseMatrixDensifyQuickTest<float> dense_matrix_densify_quick_test_float("float");
DenseMatrixDensifyQuickTest<double> dense_matrix_densify_quick_test_double("double");

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
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm0(size, size, DataType_(1));
                DenseMatrix<DataType_> dm1(size, size, DataType_(1));

                TEST_CHECK_EQUAL(dm0, dm1);
            }
        }
};
DenseMatrixEqualityTest<float> dense_matrix_equality_test_float("float");
DenseMatrixEqualityTest<double> dense_matrix_equality_test_double("double");

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
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);

                DenseMatrix<DataType_> dm(rows, columns);
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()), i_end(dm.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(i.index());
                }


                TEST_CHECK_EQUAL(dm.columns(), columns);
                TEST_CHECK_EQUAL(dm.rows(), rows);

                DenseVectorRange<DataType_> row1(dm[0]);
                DenseVector<DataType_> col1(dm.column(0));

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
            DenseMatrix<DataType_> dm(rows, columns, DataType_(1));

            TEST_CHECK_EQUAL(dm, dm);
            TEST_CHECK_EQUAL(dm.columns(), columns);
            TEST_CHECK_EQUAL(dm.rows(), rows);

            DenseVectorRange<DataType_> row1((dm)[0]);
            DenseVector<DataType_> col1(dm.column(0));

            TEST_CHECK_EQUAL(row1.size(), columns);
            TEST_CHECK_EQUAL(col1.size(), rows);

            DenseMatrix<DataType_> dm2(4, 3, DataType_(2));
            DenseMatrix<DataType_> dm3(5, 3, DataType_(2));
            DenseMatrix<DataType_> dm4(5, 4, DataType_(2));
            TEST_CHECK_THROWS(dm2 == dm3, MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(dm3 == dm4, MatrixColumnsDoNotMatch);

        }
};
DenseMatrixQuickTest<float>  dense_matrix_quick_test_float("float");
DenseMatrixQuickTest<double> dense_matrix_quick_test_double("double");
