/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Thorsten Deinert <thosten.deinert@uni-dortmund.de>
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


#include <honei/graph/evolving_graph.hh>
#include <honei/graph/position.hh>
#include <honei/graph/evolving_animator.hh>

#include <string>
#include <iostream>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
// our PositionsTest is a BaseTest.
class EvolvingGraphTest :
    public BaseTest
 {
    private:
        typedef Node<DataType_> NodeType;
    public:
        // overridden to draw a nice hello when test starts...
        EvolvingGraphTest(const std::string & type) :
            BaseTest("evolving_graph_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

    virtual void run() const
         {
            /* // create evolving graph and some nodes
            EvolvingGraph<DataType_> eg(2,1);
            eg.add_node(new NodeType(1, 1, 1, 1));
            eg.add_node(new NodeType(2, 2, 2, 2));
            eg.add_node(new NodeType(3, 3, 3, 3));

            // create a graph for timeslice 1, add some nodes already defined and define some edges
            Graph<DataType_> * t1 = new Graph<DataType_>(3,2);
            t1->add_node(eg.get_node(1));
            t1->add_node(eg.get_node(2));
            t1->add_node(eg.get_node(3));
            t1->add_edge(1,2, 1);
            t1->add_edge(2,3, 2);
            // add the whole stuff to the ev.graph
            eg.addTimeslice(t1);

            // same procedure for timeslice 2
            eg.add_node(new NodeType(4, 4, 4, 4));
            Graph<DataType_> * t2 = new Graph<DataType_>(3,2);
            t2->add_node(eg.get_node(1));
            t2->add_node(eg.get_node(3));
            t2->add_node(eg.get_node(4));
            t2->add_edge(3,4, 3);
            t2->add_edge(1,4, 4);
            eg.addTimeslice(t2);

            // and also for timeslice 3
            eg.add_node(new NodeType(5, 5, 5, 5));
            eg.add_node(new NodeType(6, 6, 6, 6));
            t2 = new Graph<DataType_>(5,2);
            t2->add_node(eg.get_node(1));
            t2->add_node(eg.get_node(3));
            t2->add_node(eg.get_node(4));
            t2->add_node(eg.get_node(5));
            t2->add_node(eg.get_node(6));
            t2->add_edge(3,4, 3);
            t2->add_edge(1,4, 4);
            t2->add_edge(3,5, 2);
            t2->add_edge(4,5, 4);
            eg.addTimeslice(t2);*/
            //std::cout << "90\n";
            EvolvingGraph<DataType_> eg(2,1);
            // Creating nodes

            //std::cout << "94\n";
            Graph<DataType_> & t1(eg.add_timeslice(2));
            Graph<DataType_> & t2(eg.add_timeslice(3));

            //std::cout << "98\n";
            t1.add_node(1, 1);
            t1.add_node(4, 1);
            //std::cout << "101\n";
            t1.add_edge(1, 4, 1);

            t2.add_node(1);
            t2.add_node(2, 2);
            t2.add_node(3, 3);
            t2.add_edge(1,2,1);
            t2.add_edge(1,3,1);

            // see what happens if we assemble the whole graph's matrices
            std::cout << "evolving graph with " << eg.node_count() << " nodes and " << eg.slice_count() << " timeslices " << std::endl;
            std::cout << "Coordinates: " << std::endl;
            std::cout << *eg.coordinates() << std::endl;

            std::cout << " NodeWeights: " << std::endl;
            std::cout << *eg.node_weights() << std::endl;
            std::cout << "Edges: " << std::endl;
            std::cout << *eg.edges() << std::endl;

            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> positions(eg, (DataType_)2);
            positions.update(0.01, 100);
            std::cout << "coordinates for eg after position.update(0.01, 100):\n" << positions.coordinates() << std::endl;

            eg.update_slice_coordinates();
            for (int i(0); i < eg.slice_count(); ++i)
                std::cout << "timeslice " << i << " after update:\n" << *eg.get_timeslice(i).coordinates() << std::endl;

            //std::string file("evolving.gml");
            //eg.write_gml(file.c_str(), true);


            EvolvingAnimator<Tag_, DataType_> animator(eg, 0.1f);

            for (animator.prepare_interpolation(); !animator.end(); animator.next_step())
                std::cout << animator.time() << std::endl;
            TEST_CHECK(true);
    }
};

// instantiate test cases
EvolvingGraphTest<tags::CPU, float> evolving_graph_test_float("float");
EvolvingGraphTest<tags::CPU, double> evolving_graph_test_doube("double");
#ifdef HONEI_SSE
EvolvingGraphTest<tags::CPU::SSE, float> sse_evolving_graph_test_float("SSE float");
EvolvingGraphTest<tags::CPU::SSE, double> sse_evolving_graph_test_doube("SSE double");
#endif
#ifdef HONEI_CELL
EvolvingGraphTest<tags::Cell, float> cell_evolving_graph_test_float("cell float");
#endif
