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

template <typename DataType_, typename Tag_>
class BreadthFirstSearchQuickTest :
    public QuickTest
{
    public:
        BreadthFirstSearchQuickTest(const std::string & type) :
            QuickTest("breadth_first_search_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
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
            unsigned long i(0);
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
            DenseMatrix<long int> distance(7, 7, 0);
            bool coherent(BreadthFirstSearch<Tag_>::value(distance, pAdj));
            //std::cout<< "distance" << distance << std::endl;

            TEST_CHECK(coherent);
            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);
            TEST_CHECK_EQUAL(distance[0][4], 3);
            TEST_CHECK_EQUAL(distance[3][4], 2);
            TEST_CHECK_EQUAL(distance[4][0], 3);
            TEST_CHECK_EQUAL(distance[4][3], 2);

            SparseMatrix<bool> pAdj2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance, pAdj2), MatrixIsNotSquare);

            SparseMatrix<bool> pAdj3(2, 2);
            pAdj3[0][1] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pAdj3[1][0] = true;
            pAdj3[1][1] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<long int> distance_2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_2, pAdj), MatrixIsNotSquare);

            DenseMatrix<long int> distance_3(6, 6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_3, pAdj), MatrixRowsDoNotMatch);

            SparseMatrix<bool>  pAdj4(7,7);
            adj[1] = 0;
            adj[7] = 0;
            i = 0;
            for (typename MutableMatrix<bool>::ElementIterator e(pAdj4.begin_elements()), e_end(pAdj4.end_elements()); e != e_end ; ++e)
            {
                if (adj[i] != 0)
                {
                    *e = adj[i];
                }
                i++;
            }

            coherent = BreadthFirstSearch<Tag_>::value(distance, pAdj4);
            TEST_CHECK(!coherent);


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
            coherent = BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW);

            TEST_CHECK(coherent);
            TEST_CHECK_EQUAL(distance2[0][0], 0);
            TEST_CHECK_EQUAL(distance2[1][1], 0);
            TEST_CHECK_EQUAL(static_cast<long int>(distance2[0][4] * 100000), 386834);
            TEST_CHECK_EQUAL(static_cast<long int>(distance2[3][4] * 100000), 366854);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance2[4][0], distance2[0][4], 2 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance2[4][3], distance2[3][4], 2 * std::numeric_limits<DataType_>::epsilon());

            SparseMatrix<DataType_>  pEW2(7,6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW2), MatrixIsNotSquare);

            SparseMatrix<DataType_> pEW3(7, 7);
            pEW3[0][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pEW3[1][0] = 4.3;
            pEW3[1][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<DataType_> distance3(7, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance3, pNW, pEW), MatrixIsNotSquare);

            DenseMatrix<DataType_> distance4(8, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance4, pNW, pEW), MatrixRowsDoNotMatch);

            DenseVector<DataType_>  pNW2(8);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW2, pEW), VectorSizeDoesNotMatch);

            DenseVector<DataType_>  pNW3(7);
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW3, pEW);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            i = 0;
            SparseMatrix<DataType_>  pEW4(7,7);
            EW[1] = 0;
            EW[7] = 0;
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEW.begin_elements()), e_end(pEW.end_elements()); e != e_end ; ++e)
            {
                if (EW[i] != 0)
                {
                    *e = 2 * EW[i];
                }
                i++;
            }
            coherent = BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW4);
            TEST_CHECK(!coherent);
        }
};

// instantiate test cases
BreadthFirstSearchQuickTest<float, tags::CPU> breadth_first_search_quick_test_float("float");
BreadthFirstSearchQuickTest<double, tags::CPU> breadth_first_search_quick_test_double("double");
BreadthFirstSearchQuickTest<float, tags::CPU::MultiCore> mc_breadth_first_search_quick_test_float("mc float");
BreadthFirstSearchQuickTest<double, tags::CPU::MultiCore> mc_breadth_first_search_quick_test_double("mc double");


template <typename DataType_, typename Tag_>
class BreadthFirstSearchCliqueTest :
    public BaseTest
{
    private:

        unsigned long _nodecount;

    public:
        BreadthFirstSearchCliqueTest(const std::string & type, unsigned long nodecount = 50) :
            BaseTest("breadth_first_search_clique_test<" + type + ">")
        {
            register_tag(Tag_::name);
            _nodecount = nodecount;
        }

        virtual void run() const
        {
            // Create a nodecount-clique
            SparseMatrix<bool>  pAdj(_nodecount, _nodecount);
            for (typename MutableMatrix<bool>::ElementIterator e(pAdj.begin_elements()), e_end(pAdj.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column())
                {
                    *e = true;
                }
            }

            // Creating a distance object with the test scenario
            DenseMatrix<long int> distance(_nodecount, _nodecount, 0);
            bool coherent(BreadthFirstSearch<Tag_>::value(distance, pAdj));

            TEST_CHECK(coherent);
            for (typename MutableMatrix<long int>::ElementIterator e(distance.begin_elements()), e_end(distance.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column()) 
                {
                    TEST_CHECK_EQUAL(*e, 1);
                }
                else TEST_CHECK_EQUAL(*e, 0);
            }

            SparseMatrix<bool> pAdj2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance, pAdj2), MatrixIsNotSquare);

            SparseMatrix<bool> pAdj3(2, 2);
            pAdj3[0][1] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pAdj3[1][0] = true;
            pAdj3[1][1] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<long int> distance_2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_2, pAdj), MatrixIsNotSquare);

            DenseMatrix<long int> distance_3(6, 6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_3, pAdj), MatrixRowsDoNotMatch);


            // Create a nodecount-clique
            SparseMatrix<DataType_>  pEW(_nodecount, _nodecount);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEW.begin_elements()), e_end(pEW.end_elements()); e != e_end ; ++e)
            {
                if (e.row() < e.column())
                {
                    *e = DataType_(2 + (e.index() % 2));
                    pEW[e.column()][e.row()] = *e;
                }
            }

            DenseVector<DataType_>  pNW(_nodecount);
            for (typename Vector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(4);
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Computing graph distance matrix and previous node matrix
            coherent = BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW);

            TEST_CHECK(coherent);
            for (typename MutableMatrix<DataType_>::ElementIterator e(distance2.begin_elements()), e_end(distance2.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column()) 
                {
                    DataType_ value(sqrt(pNW[e.row()] * pNW[e.column()]) / pEW[e.row()][e.column()]);
                    TEST_CHECK_EQUAL_WITHIN_EPS(*e, value, 2 * std::numeric_limits<DataType_>::epsilon());
                }
                else TEST_CHECK_EQUAL(*e, 0);
            }

            SparseMatrix<DataType_>  pEW2(7,6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW2), MatrixIsNotSquare);

            SparseMatrix<DataType_> pEW3(7, 7);
            pEW3[0][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pEW3[1][0] = 4.3;
            pEW3[1][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<DataType_> distance3(7, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance3, pNW, pEW), MatrixIsNotSquare);

            DenseMatrix<DataType_> distance4(8, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance4, pNW, pEW), MatrixRowsDoNotMatch);

            DenseVector<DataType_>  pNW2(8);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW2, pEW), VectorSizeDoesNotMatch);

            DenseVector<DataType_>  pNW3(_nodecount);
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW3, pEW);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

        }
};

// instantiate test cases
BreadthFirstSearchCliqueTest<float, tags::CPU> breadth_first_search_clique_test_float("clique, float", 100);
BreadthFirstSearchCliqueTest<double, tags::CPU> breadth_first_search_test_clique_double("clique, double", 100);
BreadthFirstSearchCliqueTest<float, tags::CPU::MultiCore> mc_breadth_first_search_clique_test_float("clique, mc float", 100);
BreadthFirstSearchCliqueTest<double, tags::CPU::MultiCore> mc_breadth_first_search_clique_test_double("clique, mc double", 100);


template <typename DataType_, typename Tag_>
class BreadthFirstSearchBinaryTreeTest :
    public BaseTest
{
    private:

        unsigned long _nodecount, _depth;

    public:
        BreadthFirstSearchBinaryTreeTest(const std::string & type, unsigned long depth = 10) :
            BaseTest("breadth_first_search_binary_tree_test<" + type + ">")
        {
            register_tag(Tag_::name);
            _depth = depth;
            _nodecount = (1 << (_depth + 1)) -1;
        }

        virtual void run() const
        {
            
            // Create a binary tree with deep _nodecount
            SparseMatrix<bool>  pAdj(_nodecount, _nodecount);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                pAdj[i][column] = true;
                pAdj[column][i] = true;
                pAdj[i][column + 1] = true;
                pAdj[column + 1][i] = true;
            }

            // Creating a distance object with the test scenario
            DenseMatrix<long int> distance(_nodecount, _nodecount, 0);
            bool coherent(BreadthFirstSearch<Tag_>::value(distance, pAdj));

            TEST_CHECK(coherent);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                TEST_CHECK_EQUAL(distance[i][column], 1);
                TEST_CHECK_EQUAL(distance[column][i], 1);
                TEST_CHECK_EQUAL(distance[i][column + 1], 1);
                TEST_CHECK_EQUAL(distance[column + 1][i], 1);
                TEST_CHECK_EQUAL(distance[i][i], 0);
            }
            unsigned long _nodecount_2((1 << _depth) -1);
            TEST_CHECK_EQUAL(distance[_nodecount_2][_nodecount - 1], 2* _depth);
            TEST_CHECK_EQUAL(distance[_nodecount - 1][_nodecount_2], 2* _depth);
            for (unsigned long i(1), node(1); i <= _depth; i++, node += (1 << (i-1)))
            {
                for (unsigned long j(1); j <= (1 << i); j++)
                {
                    TEST_CHECK_EQUAL(distance[0][node + j - 1], i);
                    TEST_CHECK_EQUAL(distance[node + j - 1][0], i);
                }
            }

            SparseMatrix<bool> pAdj2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance, pAdj2), MatrixIsNotSquare);

            SparseMatrix<bool> pAdj3(2, 2);
            pAdj3[0][0] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pAdj3[1][1] = true;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance, pAdj3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<long int> distance_2(3, 2);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_2, pAdj), MatrixIsNotSquare);

            DenseMatrix<long int> distance_3(6, 6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance_3, pAdj), MatrixRowsDoNotMatch);


            // Create a binary tree with deep _nodecount
            SparseMatrix<DataType_>  pEW(_nodecount, _nodecount);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                pEW[i][column] = DataType_(2);
                pEW[column][i] = DataType_(2);
                pEW[i][column + 1] = DataType_(2);
                pEW[column + 1][i] = DataType_(2);
            }

            DenseVector<DataType_>  pNW(_nodecount);
            unsigned long exp(0);
            for (typename Vector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(1 + (exp % 2));
                if (e.index() + 1 == (1 << (exp +1)) -1) exp++;
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Computing graph distance matrix and previous node matrix
            coherent = BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW);

            TEST_CHECK(coherent);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(distance2[i][column], sqrt(pNW[i] * pNW[column]) / pEW[i][column], 2 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(distance2[column][i], distance2[i][column], 2 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(distance2[i][column + 1], sqrt(pNW[i] * pNW[column + 1]) / pEW[i][column + 1], 2 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(distance2[column + 1][i], distance2[i][column + 1], 2 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL(distance[i][i], 0);
            }

            TEST_CHECK_EQUAL(static_cast<long int>(distance2[_nodecount_2][_nodecount - 1] * 100000), static_cast<long int>(2* _depth * sqrt(2) / 2 * 100000));
            TEST_CHECK_EQUAL_WITHIN_EPS(distance2[_nodecount - 1][_nodecount_2], distance2[_nodecount_2][_nodecount - 1], std::numeric_limits<DataType_>::epsilon());

            for (unsigned long i(1), node(1); i <= _depth; i++, node += (1 << (i-1)))
            {
                for (unsigned long j(1); j <= (1 << i); j++)
                {
                    TEST_CHECK_EQUAL(static_cast<long int>(distance2[0][node + j - 1] * 100000), static_cast<long int>(i * sqrt(2) / 2 *100000));
                    TEST_CHECK_EQUAL_WITHIN_EPS(distance2[node + j - 1][0], distance2[0][node + j - 1], std::numeric_limits<DataType_>::epsilon());
                }
            }

            SparseMatrix<DataType_>  pEW2(7,6);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW2), MatrixIsNotSquare);

            SparseMatrix<DataType_> pEW3(7, 7);
            pEW3[0][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            pEW3[1][0] = 4.3;
            pEW3[1][1] = 4.3;
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW3);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            DenseMatrix<DataType_> distance3(7, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance3, pNW, pEW), MatrixIsNotSquare);

            DenseMatrix<DataType_> distance4(8, 8, DataType_(0));
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance4, pNW, pEW), MatrixRowsDoNotMatch);

            DenseVector<DataType_>  pNW2(8);
            TEST_CHECK_THROWS(BreadthFirstSearch<Tag_>::value(distance2, pNW2, pEW), VectorSizeDoesNotMatch);

            DenseVector<DataType_>  pNW3(_nodecount);
            try
            {
                BreadthFirstSearch<Tag_>::value(distance2, pNW3, pEW);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

        }
};

// instantiate test cases
BreadthFirstSearchBinaryTreeTest<float, tags::CPU> breadth_first_search_binary_tree_test_float("binary_tree, float", 10);
BreadthFirstSearchBinaryTreeTest<double, tags::CPU> breadth_first_search_binary_tree_test_double("binary_tree, double", 10);
BreadthFirstSearchBinaryTreeTest<float, tags::CPU::MultiCore> mc_breadth_first_search_binary_tree_test_float("binary_tree, mc float", 10);
BreadthFirstSearchBinaryTreeTest<double, tags::CPU::MultiCore> mc_breadth_first_search_binary_tree_test_double("binary_tree, mc double", 10);
