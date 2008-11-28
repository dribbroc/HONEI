/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_matrix_tile.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using namespace tests;

template <typename DataType_>
class DenseMatrixTileCreationTest :
    public BaseTest
{
    public:
        DenseMatrixTileCreationTest(const std::string & type) :
            BaseTest("dense_matrix_tile_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm(size, size, DataType_(0));
                DenseMatrixTile<DataType_> dmt(dm, size - 2, size - 2, 2, 1);
                DenseMatrixTile<DataType_> dmt2(dmt, size - 4, size - 4, 3, 2);
                TEST_CHECK(true);
            }
        }
};

DenseMatrixTileCreationTest<float> dense_matrix_tile_creation_test_float("float");
DenseMatrixTileCreationTest<double> dense_matrix_tile_creation_test_double("double");

template <typename DataType_>
class DenseMatrixTileCopyTest :
    public BaseTest
{
    public:
        DenseMatrixTileCopyTest(const std::string & type) :
            BaseTest("dense_matrix_tile_copy_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 4) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm(size, size, DataType_(0));
                DenseMatrixTile<DataType_> dmt(dm, size - 3, size - 3, 1, 2);
                DenseMatrix<DataType_> c(dmt.copy());

                for (typename DenseMatrix<DataType_>::ElementIterator i(c.begin_elements()),
                        i_end(c.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    *i = 1;
                }

                for (typename DenseMatrixTile<DataType_>::ConstElementIterator i(dmt.begin_elements()),
                        i_end(dmt.end_elements()) ; i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
DenseMatrixTileCopyTest<float> dense_matrix_tile_copy_test_float("float");
DenseMatrixTileCopyTest<double> dense_matrix_tile_copy_test_double("double");

template <typename DataType_>
class DenseMatrixTileEqualityTest :
    public BaseTest
{
    public:
        DenseMatrixTileEqualityTest(const std::string & type) :
            BaseTest("dense_matrix_tile_equality_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(20) ; size < (1 << 7) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm0(size, size, DataType_(0));
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm0.begin_elements()), i_end(dm0.end_elements());
                        i != i_end; ++i)
                {
                    *i += i.column() + i.row();
                }

                DenseMatrixTile<DataType_> dmt0(dm0, size / 10, size / 10, 0, size / 5);
                DenseMatrixTile<DataType_> dmt1(dm0, size / 10, size / 10, size / 5, 0);
                DenseMatrixTile<DataType_> dmt2(dmt0, size / 20, size / 20, 0, size / 10);
                DenseMatrixTile<DataType_> dmt3(dmt1, size / 20 , size / 20, size / 10, 0);
                TEST_CHECK_EQUAL(dmt0, dmt1);
                TEST_CHECK_EQUAL(dmt2, dmt3);
            }
        }
};
DenseMatrixTileEqualityTest<float> dense_matrix_tile_equality_test_float("float");
DenseMatrixTileEqualityTest<double> dense_matrix_tile_equality_test_double("double");

template <typename DataType_>
class DenseMatrixTileLayoutTest :
    public BaseTest
{
    public:
        DenseMatrixTileLayoutTest(const std::string & type) :
            BaseTest("dense_matrix_tile_layout_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                unsigned long columns(size + 1), rows(size);

                DenseMatrix<DataType_> dm(rows, columns);
                DenseMatrixTile<DataType_> dmt(dm, size - 6, size - 8, 3, 7);
                for (typename DenseMatrix<DataType_>::ElementIterator i(dm.begin_elements()), i_end(dm.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = DataType_(i.index());
                }


                TEST_CHECK_EQUAL(dm.columns(), columns);
                TEST_CHECK_EQUAL(dm.rows(), rows);
                TEST_CHECK_EQUAL(dmt.columns(), size - 8);
                TEST_CHECK_EQUAL(dmt.rows(), size - 6);

                typename DenseMatrix<DataType_>::ConstRow row1(dm[0]);
                typename DenseMatrix<DataType_>::ConstColumn col1(dm.column(0));

                DenseVectorRange<DataType_> dmtrow1(dmt[0]);
                DenseVectorSlice<DataType_> dmtcol1(dmt.column(0));

                TEST_CHECK_EQUAL(row1.size(), columns);
                TEST_CHECK_EQUAL(col1.size(), rows);
                TEST_CHECK_EQUAL(dmtrow1.size(), size - 8);
                TEST_CHECK_EQUAL(dmtcol1.size(), size - 6);

                for (typename DenseMatrix<DataType_>::ConstRow::ConstElementIterator i(row1.begin_elements()), i_end(row1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, i.index(), std::numeric_limits<DataType_>::epsilon());
                }

                for (typename DenseMatrix<DataType_>::ConstColumn::ConstElementIterator i(col1.begin_elements()), i_end(col1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, i.index() * columns, std::numeric_limits<DataType_>::epsilon());
                }

                for (typename DenseMatrix<DataType_>::Row::ConstElementIterator i(dmtrow1.begin_elements()), i_end(dmtrow1.end_elements());
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, 3 * columns + 7 + i.index(), std::numeric_limits<DataType_>::epsilon());
                }

                for (typename DenseVectorSlice<DataType_>::ConstElementIterator i(dmtcol1.begin_elements()), i_end(dmtcol1.end_elements());
                        i != i_end ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, (3 + i.index()) * columns + 7, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};

DenseMatrixTileLayoutTest<float> dense_matrix_tile_layout_test_float("float");
DenseMatrixTileLayoutTest<double> dense_matrix_tile_layout_test_double("double");

template <typename DataType_>
class DenseMatrixTileQuickTest :
    public QuickTest
{
    public:
        DenseMatrixTileQuickTest(const std::string & type) :
            QuickTest("dense_matrix_tile_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long columns(3), rows(2);
            DenseMatrix<DataType_> dm(rows, columns, DataType_(1));
            DenseMatrixTile<DataType_> dmt(dm, rows - 1, columns - 1, 0, 1);

            TEST_CHECK_EQUAL(dm, dm);
            TEST_CHECK_EQUAL(dmt, dmt);
            TEST_CHECK_EQUAL(dm.columns(), columns);
            TEST_CHECK_EQUAL(dm.rows(), rows);
            TEST_CHECK_EQUAL(dmt.rows(), rows - 1);
            TEST_CHECK_EQUAL(dmt.columns(), columns - 1);

            typename DenseMatrix<DataType_>::ConstRow row1(dm[0]);
            typename DenseMatrix<DataType_>::ConstColumn col1(dm.column(0));

            DenseVectorRange<DataType_> dmtrow1(dmt[0]);
            DenseVectorSlice<DataType_> dmtcol1(dmt.column(0));

            TEST_CHECK_EQUAL(row1.size(), columns);
            TEST_CHECK_EQUAL(col1.size(), rows);
            TEST_CHECK_EQUAL(dmtrow1.size(), columns - 1);
            TEST_CHECK_EQUAL(dmtcol1.size(), rows - 1);

            dmt(1, 1)= DataType_(5);
            TEST_CHECK_EQUAL_WITHIN_EPS(dm(1, 2), DataType_(5), std::numeric_limits<DataType_>::epsilon());

            DenseMatrix<DataType_> dm2(4, 3, DataType_(2));
            DenseMatrix<DataType_> dm3(5, 3, DataType_(2));
            DenseMatrixTile<DataType_> dmt2(dm2, 2, 2, 0, 0);
            DenseMatrixTile<DataType_> dmt3(dm3, 1, 1, 1, 1);
            TEST_CHECK_THROWS(dmt2==dmt, MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(dmt3==dmt, MatrixColumnsDoNotMatch);

            DenseMatrix<DataType_> dm4(7, 7, DataType_(19));
            DenseMatrixTile<DataType_> dmt4_1(dm4, 4, 4, 0, 3);
            DenseMatrixTile<DataType_> dmt4_2(dm4, 4, 4, 3, 0);

            for (typename DenseMatrixTile<DataType_>::ElementIterator i(dmt4_1.begin_elements()), i_end(dmt4_1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(26);
            }

            TEST_CHECK_EQUAL(dmt4_1(3, 0), DataType_(26));
            dmt4_1(3, 0) = DataType_(19);

            for (typename DenseMatrixTile<DataType_>::ConstElementIterator i(dmt4_2.begin_elements()), i_end(dmt4_2.end_elements());
                    i != i_end; ++i)
            {
                TEST_CHECK_EQUAL(*i, DataType_(19));
            }

        }
};
DenseMatrixTileQuickTest<float>  dense_matrix_tile_quick_test_float("float");
DenseMatrixTileQuickTest<double> dense_matrix_tile_quick_test_double("double");
