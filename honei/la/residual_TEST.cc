/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/residual.hh>
#include <honei/la/product.hh>
#include <honei/util/unittest.hh>

#include <limits>

using namespace honei;
using namespace tests;

template <typename DataType_>
class DenseResidualTest :
    public BaseTest
{
    public:
        DenseResidualTest(const std::string & type) :
            BaseTest("dense_residual_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                DataType_ s(size);

                {
                    DenseVector<DataType_> b(size), x(size);
                    DenseMatrix<DataType_> a(size, size);

                    for (typename DenseVector<DataType_>::ElementIterator i(x.begin_elements()),
                            i_end(x.end_elements()), j(b.begin_elements()) ; i != i_end ; ++i, ++j)
                    {
                        *i = DataType_(i.index() + 1);
                        *j = DataType_(j.index()) * s * (s + 1) / 2;
                    }

                    for (typename DenseMatrix<DataType_>::ElementIterator i(a.begin_elements()),
                            i_end(a.end_elements()) ; i != i_end ; ++i)
                    {
                        *i = DataType_(i.row());
                    }

                    Residual<>::value(b, a, x);

                    for (typename DenseVector<DataType_>::ConstElementIterator i(b.begin_elements()),
                            i_end(b.end_elements()) ; i != i_end ; ++i)
                    {
                        DataType_ r(i.index() * s * s * s); /// \todo Find a lower border.
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, r * std::numeric_limits<DataType_>::epsilon());
                    }
                }

                DenseVector<DataType_> b1(size + 1), b2(size);
                DenseVector<DataType_> x1(size), x2(size + 1);
                DenseMatrix<DataType_> a1(size + 1, size + 1), a2(size+1, size);
                TEST_CHECK_THROWS(Residual<>::value(b1, a1, x1), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b2, a1, x2), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b1, a2, x2), MatrixIsNotSquare);
            }
        }
};
DenseResidualTest<float> dense_residual_test_float("float");
DenseResidualTest<double> dense_residual_test_double("double");

template <typename DataType_>
class DenseResidualQuickTest :
    public QuickTest
{
    public:
        DenseResidualQuickTest(const std::string & type) :
            QuickTest("dense_residual_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DataType_ s(size);

            {
                DenseVector<DataType_> b(size), x(size);
                DenseMatrix<DataType_> a(size, size);

                for (typename DenseVector<DataType_>::ElementIterator i(x.begin_elements()),
                        i_end(x.end_elements()), j(b.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    *i = DataType_(i.index() + 1);
                    *j = DataType_(j.index()) * s * (s + 1) / 2;
                }

                for (typename DenseMatrix<DataType_>::ElementIterator i(a.begin_elements()),
                        i_end(a.end_elements()) ; i != i_end ; ++i)
                {
                    *i = DataType_(i.row());
                }

                Residual<>::value(b, a, x);

                for (typename DenseVector<DataType_>::ConstElementIterator i(b.begin_elements()),
                        i_end(b.end_elements()) ; i != i_end ; ++i)
                {
                    DataType_ r(i.index() * s * s * s); /// \todo Find a lower border.
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, r * std::numeric_limits<DataType_>::epsilon());
                }

                DenseVector<DataType_> b1(size + 1), b2(size);
                DenseVector<DataType_> x1(size), x2(size + 1);
                DenseMatrix<DataType_> a1(size + 1, size + 1), a2(size+1, size);
                TEST_CHECK_THROWS(Residual<>::value(b1, a1, x1), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b2, a1, x2), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b1, a2, x2), MatrixIsNotSquare);
            }
        }
};
DenseResidualQuickTest<float> dense_residual_quick_test_float("float");
DenseResidualQuickTest<double> dense_residual_quick_test_double("double");


template <typename DataType_>
class SparseResidualTest :
    public BaseTest
{
    public:
        SparseResidualTest(const std::string & type) :
            BaseTest("sparse_residual_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 6) ; size <<= 1)
            {
                DataType_ s(size);

                {
                    DenseVector<DataType_> b(size), x(size);
                    SparseMatrix<DataType_> a(size, size);

                    for (typename DenseVector<DataType_>::ElementIterator i(x.begin_elements()),
                            i_end(x.end_elements()), j(b.begin_elements()) ; i != i_end ; ++i, ++j)
                    {
                        *i = DataType_(i.index() + 1);
                        *j = DataType_(j.index()) * s * (s + 1) / 2;
                    }

                    for (typename SparseMatrix<DataType_>::ElementIterator i(a.begin_elements()),
                            i_end(a.end_elements()) ; i != i_end ; ++i)
                    {
                        *i = DataType_(i.row());
                    }

                    Residual<>::value(b, a, x);

                    for (typename DenseVector<DataType_>::ConstElementIterator i(b.begin_elements()),
                            i_end(b.end_elements()) ; i != i_end ; ++i)
                    {
                        DataType_ r(i.index() * s * s * s); /// \todo Find a lower border.
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, r * std::numeric_limits<DataType_>::epsilon());
                    }
                }

                DenseVector<DataType_> b1(size + 1), b2(size);
                DenseVector<DataType_> x1(size), x2(size + 1);
                SparseMatrix<DataType_> a1(size + 1, size + 1), a2(size+1, size);
                TEST_CHECK_THROWS(Residual<>::value(b1, a1, x1), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b2, a1, x2), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b1, a2, x2), MatrixIsNotSquare);
            }
        }
};
SparseResidualTest<float>  sparse_residual_test_float("float");
SparseResidualTest<double> sparse_residual_test_double("double");


template <typename DataType_>
class SparseResidualQuickTest :
    public QuickTest
{
    public:
        SparseResidualQuickTest(const std::string & type) :
            QuickTest("sparse_residual_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DataType_ s(size);

            {
                DenseVector<DataType_> b(size), x(size);
                SparseMatrix<DataType_> a(size, size);

                for (typename DenseVector<DataType_>::ElementIterator i(x.begin_elements()),
                        i_end(x.end_elements()), j(b.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    *i = DataType_(i.index() + 1);
                    *j = DataType_(j.index()) * s * (s + 1) / 2;
                }

                for (typename SparseMatrix<DataType_>::ElementIterator i(a.begin_elements()),
                        i_end(a.end_elements()) ; i != i_end ; ++i)
                {
                    *i = DataType_(i.row());
                }

                Residual<>::value(b, a, x);

                for (typename DenseVector<DataType_>::ConstElementIterator i(b.begin_elements()),
                        i_end(b.end_elements()) ; i != i_end ; ++i)
                {
                    DataType_ r(i.index() * s * s * s); /// \todo Find a lower border.
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, r * std::numeric_limits<DataType_>::epsilon());
                }

                DenseVector<DataType_> b1(size + 1), b2(size);
                DenseVector<DataType_> x1(size), x2(size + 1);
                DenseMatrix<DataType_> a1(size + 1, size + 1), a2(size+1, size);
                TEST_CHECK_THROWS(Residual<>::value(b1, a1, x1), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b2, a1, x2), VectorSizeDoesNotMatch);
                TEST_CHECK_THROWS(Residual<>::value(b1, a2, x2), MatrixIsNotSquare);
            }
        }
};
SparseResidualQuickTest<float>  sparse_residual_quick_test_float("float");
SparseResidualQuickTest<double> sparse_residual_quick_test_double("double");
