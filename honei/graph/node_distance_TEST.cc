/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 *
 * This file is part of the Graph C++ library. LibGraph is free software;
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

#include <honei/la/dense_matrix.hh>
#include <honei/graph/node_distance.hh>
#include <honei/util/unittest.hh>

#include <string>
#include <iostream>
#include <math.h>

using namespace honei;
using  namespace tests;
template <typename DataType_,typename Tag_>
class NodeDistanceQuickTest :
    public QuickTest
{
    public:
        NodeDistanceQuickTest(const std::string & type) :
            QuickTest("node_distance_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long columns(2), rows(3);
            DenseMatrix<DataType_> dm(rows, columns);

            for (typename DenseMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                *i = i.index();
            }


            DenseMatrix<DataType_> weight_of_edges(3, 3, DataType_(0));
            weight_of_edges[0][1] = 1;
            weight_of_edges[1][0] = 4;
            weight_of_edges[1][2] = 3;
            weight_of_edges[2][1] = 2;

            DenseMatrix<DataType_> pSquarDist(3, 3, DataType_(0));
            DenseMatrix<DataType_> pInv_SquarDist(3, 3, DataType_(0));
            DataType_ range(3);


            DenseMatrix<DataType_> distance = NodeDistance<Tag_>::value(dm);

            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][0], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][2], 32, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][2], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][0], 32, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][2], 0, std::numeric_limits<DataType_>::epsilon());


            NodeDistance<Tag_>::value(dm, weight_of_edges, pSquarDist, pInv_SquarDist, range);

            TEST_CHECK_EQUAL((pSquarDist)[0][0], 0);
            TEST_CHECK_EQUAL((pSquarDist)[1][1], 0);
            TEST_CHECK_EQUAL((pSquarDist)[2][2], 0);
            TEST_CHECK_EQUAL((pSquarDist)[2][0], 0);
            TEST_CHECK_EQUAL((pSquarDist)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((pSquarDist)[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pSquarDist)[1][0], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pSquarDist)[1][2], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pSquarDist)[2][1], 8, std::numeric_limits<DataType_>::epsilon());

            TEST_CHECK_EQUAL((pInv_SquarDist)[0][0], 0);
            TEST_CHECK_EQUAL((pInv_SquarDist)[1][1], 0);
            TEST_CHECK_EQUAL((pInv_SquarDist)[2][2], 0);
            TEST_CHECK_EQUAL((pInv_SquarDist)[2][0], 0);
            TEST_CHECK_EQUAL((pInv_SquarDist)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((pInv_SquarDist)[0][1], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pInv_SquarDist)[1][0], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pInv_SquarDist)[1][2], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((pInv_SquarDist)[2][1], 0.125, std::numeric_limits<DataType_>::epsilon());
        }
};
NodeDistanceQuickTest<float, tags::CPU>  node_distance_quick_test_float("float");
NodeDistanceQuickTest<double, tags::CPU> node_distance_quick_test_double("double");
NodeDistanceQuickTest<float, tags::CPU::MultiCore>  mc_node_distance_quick_test_float("mc float");
NodeDistanceQuickTest<double, tags::CPU::MultiCore>  mc_node_distance_quick_test_double("mc double");
#ifdef HONEI_SSE
NodeDistanceQuickTest<float, tags::CPU::SSE>  sse_node_distance_quick_test_float("sse float");
NodeDistanceQuickTest<double, tags::CPU::SSE> sse_node_distance_quick_test_double("sse double");
NodeDistanceQuickTest<float, tags::CPU::MultiCore::SSE>  mc_sse_node_distance_quick_test_float("mc sse float");
NodeDistanceQuickTest<double, tags::CPU::MultiCore::SSE> mc_sse_node_distance_quick_test_double("mc sse double");
#endif
#ifdef HONEI_CELL
NodeDistanceQuickTest<float, tags::Cell> cell_node_distance_quick_test_float("Cell float");
#endif

template <typename DataType_,typename Tag_>
class NodeDistanceTest :
    public BaseTest
{
    private:
        unsigned long _nodecount;
    public:
        NodeDistanceTest(const std::string & type, unsigned long nodecount) :
            BaseTest("node_distance_test<" + type + ">")
        {
            register_tag(Tag_::name);
            _nodecount = nodecount;
        }

        virtual void run() const
        {
            unsigned long columns(2), rows(_nodecount);
            DenseMatrix<DataType_> dm(rows, columns, DataType_(0));

            for (typename DenseMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                i.index() % 2 == 0 ? *i = i.index() : *i = 0;
            }

            DenseMatrix<DataType_> weight_of_edges(_nodecount, _nodecount, DataType_(0));
            for (typename DenseMatrix<DataType_>::ElementIterator i(weight_of_edges.begin_elements()),
                    i_end(weight_of_edges.end_elements()) ; i != i_end ; ++i)
            {
                i.index()%3 != 0 ? *i = i.index() / 13 : 0;
            }

            DenseMatrix<DataType_> pSquarDist(_nodecount, _nodecount, DataType_(0));
            DenseMatrix<DataType_> pInv_SquarDist(_nodecount, _nodecount, DataType_(0));
            DataType_ range(1.5);

            DenseMatrix<DataType_> distance = NodeDistance<Tag_>::value(dm);

            for (typename DenseMatrix<DataType_>::ElementIterator i(distance.begin_elements()),
                    i_end(distance.end_elements()) ; i != i_end ; ++i)
            {
                DataType_ value(abs(i.row()-i.column()) * 2);
                value *= value;
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, value, std::numeric_limits<DataType_>::epsilon());
            }

            NodeDistance<Tag_>::value(dm, weight_of_edges, pSquarDist, pInv_SquarDist, range);

            typename DenseMatrix<DataType_>::ElementIterator j(pSquarDist.begin_elements());
            typename DenseMatrix<DataType_>::ElementIterator g(weight_of_edges.begin_elements());
            for (typename DenseMatrix<DataType_>::ElementIterator i(pInv_SquarDist.begin_elements()),
                    i_end(pInv_SquarDist.end_elements()) ; i != i_end ; ++i, ++j, ++g)
            {
                DataType_ value(abs(i.row()-i.column()) * 2);
                value *= value;
                DataType_ value_1(0);
                DataType_ value_2(0);
                (value < range) && (value > std::numeric_limits<DataType_>::epsilon()) ? value_1 = 1 / value : value_1 = 0;
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, value_1, std::numeric_limits<DataType_>::epsilon());
                if (*g > std::numeric_limits<DataType_>::epsilon()) value_2 = value;
                TEST_CHECK_EQUAL_WITHIN_EPS(*j, value_2, std::numeric_limits<DataType_>::epsilon());
            }
        }
};
NodeDistanceTest<float, tags::CPU> node_distance_test_float("float", 10);
NodeDistanceTest<double, tags::CPU> node_distance_test_double("double", 10);
NodeDistanceTest<float, tags::CPU::MultiCore>  mc_node_distance_test_float("mc float", 10);
NodeDistanceTest<double, tags::CPU::MultiCore>  mc_node_distance_test_double("mc double", 10);
#ifdef HONEI_SSE
NodeDistanceTest<float, tags::CPU::SSE>  sse_node_distance_test_float("sse float", 10);
NodeDistanceTest<double, tags::CPU::SSE> sse_node_distance_test_double("sse double", 10);
NodeDistanceTest<float, tags::CPU::MultiCore::SSE>  mc_sse_node_distance_test_float("mc sse float", 10);
NodeDistanceTest<double, tags::CPU::MultiCore::SSE> mc_sse_node_distance_test_double("mc sse double", 10);
#endif
#ifdef HONEI_CELL
NodeDistanceTest<float, tags::Cell> cell_node_distance_test_float("Cell float", 10);
#endif
