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

#include <honei/util/unittest.hh>


#include <honei/graph/graph.hh>
#include <honei/graph/abstract_graph.hh>
#include <honei/graph/position.hh>
#include <honei/graph/test_scenario.hh>

#include <string>
#include <iostream>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
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
        Graph<DataType_> * g(TestScenario<DataType_>::Grid(6, 6));
        std::cout << *g->coordinates() << std::endl;
        std::cout << *g->node_weights() << std::endl;
        std::cout << *g->edges() << std::endl;

        Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> p(*g, DataType_(1));
        p.update(0, 10000);
        //std::string file("test.gml");
        //g->write_gml(file.c_str(), true);
        TEST_CHECK(true);
        delete g;
    }
};

// instantiate test cases
#ifdef HONEI_SSE
GraphTest<tags::CPU::SSE, float> graph_test_float("sse float");
GraphTest<tags::CPU::SSE, double> graph_test_doube("sse double");
#endif
#ifdef HONEI_CELL
GraphTest<tags::Cell, float> graph_test_float("cell float");
#endif
