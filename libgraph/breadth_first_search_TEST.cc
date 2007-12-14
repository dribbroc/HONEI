/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <libgraph/breadth_first_search.hh>
#include <unittest/unittest.hh>
#include <libla/dense_matrix.hh>

#include <string>

using namespace honei;
using  namespace tests;

template <typename DataType_>
class BreadthFirstSearchQuickTest :
    public QuickTest
{
    public:
        BreadthFirstSearchQuickTest(const std::string & type) :
            QuickTest("breadth_first_search_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            bool adj[] =  {0, 1, 0, 0, 0, 0, 0,
                           1, 0, 1, 0, 0, 0, 0,
                           0, 1, 0, 1, 1, 0, 0,
                           0, 0, 1, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 1, 0,
                           0, 0, 0, 0, 1, 0, 1,
                           0, 0, 0, 0, 0, 1, 0};

            // Now, fill that numbers into the real matrices
            int i(0);
            SparseMatrix<bool>  pAdj(7,7);
            for (typename MutableMatrix<bool>::ElementIterator e(pAdj.begin_elements()), e_end(pAdj.end_elements()); e != e_end ; ++e)
                    {
                        if (adj[i] != 0)
                        {
                            *e = adj[i];
                        }
                        i++;
                    }

            // Creating a distance object with the test scenario
            DenseMatrix<int> distance(7, 7, int(0));
            bool connected(BreadthFirstSearch<tags::CPU>::value(distance, pAdj));


            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][4], 3, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[3][4], 2, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[4][0], 3, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[4][3], 2, std::numeric_limits<DataType_>::epsilon());

            SparseMatrix<bool> pAdj2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<tags::CPU>::value(distance, pAdj2), MatrixRowsDoNotMatch);


            // Creating test scenario
            DataType_ EW[]=  {0, 1, 0, 0, 0, 0, 0,
                              1, 0, 1, 0, 0, 0, 0,
                              0, 1, 0, 1, 1, 0, 0,
                              0, 0, 1, 0, 0, 0, 0,
                              0, 0, 1, 0, 0, 1, 0,
                              0, 0, 0, 0, 1, 0, 1,
                              0, 0, 0, 0, 0, 1, 0};

            DataType_ NW[]=  {1, 2, 3, 4, 5, 6, 7};

            // Now, fill that numbers into the real matrices
            i = 0;
            SparseMatrix<DataType_>  pEW(7,7);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEW.begin_elements()), e_end(pEW.end_elements()); e != e_end ; ++e)
                    {
                        if (EW[i] != 0)
                        {
                            *e = 2 * EW[i];
                        }
                        i++;
                    }

            i = 0;
            DenseVector<DataType_>  pNW(7);
            for (typename Vector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
                    {
                        *e = NW[i];
                        i++;
                    }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(7, 7, DataType_(0));

            // Computing graph distance matrix and previous node matrix
            connected = BreadthFirstSearch<tags::CPU>::value(distance2, pNW, pEW);

            TEST_CHECK_EQUAL(distance2[0][0], 0);
            TEST_CHECK_EQUAL(distance2[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance2[0][4], distance2[4][0], 2 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance2[3][4], distance2[4][3], 2 * std::numeric_limits<DataType_>::epsilon());

            DenseMatrix<DataType_> distance3(3, 7, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<tags::CPU>::value(distance3, pNW, pEW), MatrixRowsDoNotMatch);

            DenseMatrix<DataType_> distance4(7, 7, DataType_(0));
            SparseMatrix<DataType_>  pEW2(7,6);
            TEST_CHECK_THROWS(BreadthFirstSearch<tags::CPU>::value(distance4, pNW, pEW2), MatrixColumnsDoNotMatch);

            SparseMatrix<DataType_>  pEW3(6,7);
            TEST_CHECK_THROWS(BreadthFirstSearch<tags::CPU>::value(distance4, pNW, pEW3), MatrixRowsDoNotMatch);
        }
};

// instantiate test cases
BreadthFirstSearchQuickTest<float> breadth_first_search_quick_test_float("float");
BreadthFirstSearchQuickTest<double> breadth_first_search_quick_test_double("double");
