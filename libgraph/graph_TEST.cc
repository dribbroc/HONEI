/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
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

#include <unittest/unittest.hh>


#include <libgraph/graph.hh>
#include <libgraph/abstract_graph.hh>

#include <string>
#include <iostream>

using namespace honei;
using namespace tests;

template <typename DataType_>
// our PositionsTest is a BaseTest.
class GraphTest :
    public BaseTest
{
    public:
        // overridden to draw a nice hello when test starts...
        GraphTest(const std::string & type) :
            BaseTest("graph_test<" + type + ">")
        {
        }

    virtual void run() const
    {
        Graph<DataType_> g(5, 3);
        Node<DataType_> * node = new Node<DataType_>(1, 1, 1, 1, 1);
        g.addNode(node);
            g.addNode(new Node<DataType_>(2,2,2,2,2));
        g.addNode(new Node<DataType_>(3,3,3,3,3));
            g.addNode(new Node<DataType_>(4,4,4,4,4));
        g.addNode(new Node<DataType_>(6,6,6,6,6));
            g.addEdge(1,2, 7);
        g.addEdge(2,3, 4);
        g.addEdge(g.getNode(4), g.getNodeByID(6), 11);
        std::cout << *g.coordinates() << std::endl;
        std::cout << *g.nodeWeights() << std::endl;
        std::cout << *g.edges() << std::endl;

        delete(node);
        TEST_CHECK(true);
    }
};

// instantiate test cases
GraphTest<float> graph_test_float("float");
GraphTest<double> graph_test_doube("double");
