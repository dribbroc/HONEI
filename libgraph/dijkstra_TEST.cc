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


#include <libgraph/dijkstra.hh>
#include <unittest/unittest.hh>
#include <libla/dense_matrix.hh>

#include <string>

using namespace honei;
using  namespace tests;

template <typename DataType_>
class DijkstraQuickTest :
    public QuickTest
{
    public:
        DijkstraQuickTest(const std::string & type) :
            QuickTest("dijkstra_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ cost[] =  {14, 1, 14, 14, 14, 14, 14,
                                  1, 14, 1, 14, 14, 14, 14,
                                 14, 1, 14, 1, 1, 14, 14,
                                 14, 14, 1, 14, 14, 14, 14,
                                 14, 14, 1, 14, 14, 1, 14,
                                 14, 14, 14, 14, 1, 14, 1,
                                 14, 14, 14, 14, 14, 1, 14};

            // Now, fill that numbers into the real matrices
            int i(0);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pCost(new DenseMatrix<DataType_>(7,7));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pCost2(new DenseMatrix<DataType_>(7,7));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pCost->begin_elements()), f(pCost2->begin_elements()),
                    e_end(pCost->end_elements()); e != e_end ; ++e, ++f)
                    {
                        *e = cost[i];
                        *f = cost[i];
                        i++;
                    }
            
            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance = Dijkstra<DataType_>::value(*pCost);

            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][4], 3, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[3][4], 2, std::numeric_limits<DataType_>::epsilon());

            std::tr1::shared_ptr<DenseMatrix<DataType_> > pCost3(new DenseMatrix<DataType_>(3, 2));
            TEST_CHECK_THROWS(Dijkstra<DataType_>::value(*pCost3), MatrixRowsDoNotMatch);

            // Creating an empty node matrix to compute the previous node matrix
            std::tr1::shared_ptr<DenseMatrix<int> > pNodes(new DenseMatrix<int>(7,7));

            // Computing graph distance matrix and previous node matrix
            Dijkstra<DataType_>::value( *pCost2, *pNodes);

            TEST_CHECK_EQUAL((*pCost2)[0][0], 0);
            TEST_CHECK_EQUAL((*pCost2)[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((*pCost2)[0][4], 3, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pCost2)[3][4], 2, std::numeric_limits<DataType_>::epsilon());

            TEST_CHECK_EQUAL((*pNodes)[0][0], 0);
            TEST_CHECK_EQUAL((*pNodes)[2][3], 3);
            TEST_CHECK_EQUAL((*pNodes)[3][4], 2);
            TEST_CHECK_EQUAL((*pNodes)[6][4], 5);

            std::tr1::shared_ptr<DenseMatrix<DataType_> > pCost4(new DenseMatrix<DataType_>(3, 7));
            TEST_CHECK_THROWS(Dijkstra<DataType_>::value(*pCost4, *pNodes), MatrixRowsDoNotMatch);

            std::tr1::shared_ptr<DenseMatrix<int> > pNodes2(new DenseMatrix<int>(7, 3));
            TEST_CHECK_THROWS(Dijkstra<DataType_>::value(*pCost2, *pNodes2), MatrixColumnsDoNotMatch);

            std::tr1::shared_ptr<DenseMatrix<int> > pNodes3(new DenseMatrix<int>(3, 7));
            TEST_CHECK_THROWS(Dijkstra<DataType_>::value(*pCost2, *pNodes3), MatrixRowsDoNotMatch);
            
        }
};

// instantiate test cases
DijkstraQuickTest<float> dijkstra_quick_test_float("float");
DijkstraQuickTest<double> dijkstra_quick_test_double("double");
