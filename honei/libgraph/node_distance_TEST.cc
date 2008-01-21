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

#include <honei/libla/dense_matrix.hh>
#include <honei/libgraph/node_distance.hh>
#include <unittest/unittest.hh>

#include <string>
#include <iostream>

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
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(rows, columns));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm2(new DenseMatrix<DataType_>(rows, columns));

            typename MutableMatrix<DataType_>::ElementIterator e(dm2->begin_elements());

            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()),
                    i_end(dm->end_elements()) ; i != i_end ; ++i, ++e)
            {
                *i = i.index();
                *e = e.index();
            }

            std::tr1::shared_ptr<DenseMatrix<bool> > pNeighbour(new DenseMatrix<bool>(3,3, false));
            (*pNeighbour)[0][1] = true;
            (*pNeighbour)[1][0] = true;
            (*pNeighbour)[1][2] = true;
            (*pNeighbour)[2][1] = true;

            std::tr1::shared_ptr<DenseMatrix<DataType_> > weight_of_edges(new DenseMatrix<DataType_>(3,3, DataType_(0)));
            (*weight_of_edges)[0][1] = 1;
            (*weight_of_edges)[1][0] = 4;
            (*weight_of_edges)[1][2] = 3;
            (*weight_of_edges)[2][1] = 2;

            std::tr1::shared_ptr<DenseMatrix<DataType_> > pSquarDist(new DenseMatrix<DataType_>(3,3, DataType_(0)));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pInv_SquarDist(new DenseMatrix<DataType_>(3,3, DataType_(0)));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pSquarDist2(new DenseMatrix<DataType_>(3,3, DataType_(0)));
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pInv_SquarDist2(new DenseMatrix<DataType_>(3,3, DataType_(0)));
            DataType_ range(3);


            DenseMatrix<DataType_> distance = NodeDistance<Tag_>::value(*dm);

            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][0], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][2], 32, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][2], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][0], 32, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[2][2], 0, std::numeric_limits<DataType_>::epsilon());

            NodeDistance<Tag_>::value(*dm2, *pNeighbour, *pSquarDist, *pInv_SquarDist, range);

            TEST_CHECK_EQUAL((*pSquarDist)[0][0], 0);
            TEST_CHECK_EQUAL((*pSquarDist)[1][1], 0);
            TEST_CHECK_EQUAL((*pSquarDist)[2][2], 0);
            TEST_CHECK_EQUAL((*pSquarDist)[2][0], 0);
            TEST_CHECK_EQUAL((*pSquarDist)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist)[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist)[1][0], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist)[1][2], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist)[2][1], 8, std::numeric_limits<DataType_>::epsilon());

            TEST_CHECK_EQUAL((*pInv_SquarDist)[0][0], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist)[1][1], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist)[2][2], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist)[2][0], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist)[0][1], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist)[1][0], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist)[1][2], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist)[2][1], 0.125, std::numeric_limits<DataType_>::epsilon());

            NodeDistance<Tag_>::value(*dm2, *weight_of_edges, *pSquarDist2, *pInv_SquarDist2, range);

            TEST_CHECK_EQUAL((*pSquarDist2)[0][0], 0);
            TEST_CHECK_EQUAL((*pSquarDist2)[1][1], 0);
            TEST_CHECK_EQUAL((*pSquarDist2)[2][2], 0);
            TEST_CHECK_EQUAL((*pSquarDist2)[2][0], 0);
            TEST_CHECK_EQUAL((*pSquarDist2)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist2)[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist2)[1][0], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist2)[1][2], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pSquarDist2)[2][1], 8, std::numeric_limits<DataType_>::epsilon());

            TEST_CHECK_EQUAL((*pInv_SquarDist2)[0][0], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist2)[1][1], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist2)[2][2], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist2)[2][0], 0);
            TEST_CHECK_EQUAL((*pInv_SquarDist2)[0][2], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist2)[0][1], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist2)[1][0], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist2)[1][2], 0.125, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS((*pInv_SquarDist2)[2][1], 0.125, std::numeric_limits<DataType_>::epsilon());
        }
};
NodeDistanceQuickTest<float, tags::CPU>  node_distance_quick_test_float("float");
NodeDistanceQuickTest<double, tags::CPU> node_distance_quick_test_double("double");
#ifdef HONEI_SSE
NodeDistanceQuickTest<float, tags::CPU::SSE>  sse_node_distance_quick_test_float("sse float");
NodeDistanceQuickTest<double, tags::CPU::SSE> sse_node_distance_quick_test_double("sse double");
#endif
#ifdef HONEI_CELL
NodeDistanceQuickTest<float, tags::Cell> cell_node_distance_quick_test_float("Cell float");
#endif
