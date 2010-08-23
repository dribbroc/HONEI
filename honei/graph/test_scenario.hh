/* vim: set sw=4 sts=4 et nofoldenable syntax on */

/*
 * Copyright (c) 2007 Thorsten Deinert <thosten.deinert@uni-dortmund.de>
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

#pragma once
#ifndef LIBGRAPH_GUARD_TEST_SCENARIO
#define LIBGRAPH_GUARD_TEST_SCENARIO 1

#include <tr1/memory>
#include <honei/graph/graph.hh>
#include <honei/graph/evolving_graph.hh>
#include <iostream>
namespace honei
{
    template <typename DataType_> class TestScenario
    {
        public:
        static EvolvingGraph<DataType_> * Evolving(int slices, int sizes[])
        {
            EvolvingGraph<DataType_> * eg = new EvolvingGraph<DataType_>(2,  DataType_(0.5));
            for  (int t(0) ; t < slices ; ++t)
            {
                int size = sizes[t];
                int nodes(size * size);
                 std::cout << "slice " << t << ", size = " << nodes << "Nodes\n";
               
                Graph<DataType_> g(eg->add_timeslice(nodes));
                for (int n(0);  n < nodes; ++n)
                    g.add_node(n,1);
            
                for (int i = 0; i < size - 1; ++i)
                    for (int j = 0; j < size - 1; j++)
                    {
                        g.add_edge(i*size + j, i*size + j + 1);
                        g.add_edge(i*size + j, (i+1)* size + j);
                }
            for (int i = 0; i < size - 1; ++i)
                g.add_edge(i*size + size-1,  (i+1) * size + size -1);
            for (int j = 0; j < size - 1; ++j)
                g.add_edge((size-1)*size + j,  (size-1) * size + j+1);           
            }
            return eg;
        }
        
        
        static Graph<DataType_> * Grid(int rows, int columns)
        {
            int nodes = rows * columns;
            Graph<DataType_> * g = new Graph<DataType_>(nodes);
            g->random_positions(true);
            
            for (int n = 0;  n < nodes; ++n)
                g->add_node(n,1);
            
            for (int i = 0; i < rows - 1; ++i)
                for (int j = 0; j < columns - 1; j++)
                {
                    g->add_edge(i*columns + j, i*columns + j + 1);
                    g->add_edge(i*columns + j, (i+1)* columns + j);
                }
            for (int i = 0; i < rows - 1; ++i)
                g->add_edge(i*columns + columns-1,  (i+1) * columns + columns -1);
            for (int j = 0; j < columns - 1; ++j)
                g->add_edge((rows-1)*columns + j,  (rows-1) * columns + j+1);
            return g; 
        }
        
        static Graph<DataType_> * BinaryTree(int levels)
        {
            int nodes = (1 << (levels+1)) - 1;
            
            Graph<DataType_> *g = new Graph<DataType_>(nodes);
            int l = levels+1;
            int m = 2;   
            for (int n = 0; n < nodes; ++n)
            {
                g->add_node(n, l);
                if (n == m - 2)
                {
                    m *=2;
                    --l;
                }
            }
            std::cout << *g->node_weights();
            for (int n = 0; n < nodes - (1 << levels); ++n)
            {
                    g->add_edge(n, 2*n + 1);
                    g->add_edge(n, 2*n + 2);
            }
            return g;
        }
    };
}
#endif
