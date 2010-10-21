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
#ifndef LIBGRAPH_GUARD_GRAPH
#define LIBGRAPH_GUARD_GRAPH 1


#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/graph/node.hh>
#include <honei/graph/abstract_graph.hh>
#include <map>
#include <cstdlib>



namespace honei
{

    template <typename DataType_> class Graph: 
        public AbstractGraph<DataType_>
    {
    private:
        typedef    Node<DataType_> NodeType;
        NodeType** _nodes;
        int _node_count;
        int _max_nodes;
        int _coordinate_dimensions;
        std::map<int, int> _node_mapping;
        bool _random;

        /// sets edges[v_index][w_index] = edges[w_index][v_index] = weight
        void set_edge_weight_internal(int v_index, int w_index, DataType_ weight)
        {
            (*this->_edges)(v_index, w_index) = weight;
            (*this->_edges)(w_index, v_index) = weight;
        }

    public:
        /// constructs new graph with given number of nodes and a given dimension for the coordinates
        Graph(int nodes, int coordinate_dimensions=2) :
            _node_count(0),
            _max_nodes(nodes),
            _coordinate_dimensions(coordinate_dimensions),
            _node_mapping(),
            _random(false)
        {
            _nodes = new NodeType*[nodes];
            this->_coordinates = new DenseMatrix<DataType_>(nodes, coordinate_dimensions);
            this->_edges = new SparseMatrix<DataType_>(nodes, nodes);
            this->_node_weights = new DenseVector<DataType_>(nodes, (DataType_)1);
        }

        /// adds a node to this graph, puts its initial position and its weight into the relevant matrices
        void add_node(NodeType * node)
        {
            if (_node_count < _max_nodes)
            {
                _nodes[_node_count] = node;
                _node_mapping[node->id()] = _node_count;               
                DataType_ angle((DataType_)std::rand() / RAND_MAX * 3.14 * 2);
                DataType_ dist(DataType_(_max_nodes) / 20);
                if (_random)
                {
                    dist *= (DataType_)std::rand() / RAND_MAX;
                }
                (*this->_coordinates)(_node_count, 0) =  dist * sin(angle); 
                (*this->_coordinates)(_node_count, 1) =  dist * cos(angle);
                (*this->_node_weights)[_node_count] = node->get_weight();
                ++_node_count;
            }
        }
        
        inline bool random_positions(bool value)
        {
            _random = value;
	    return _random;
        }
        
        void add_node(int id, DataType_ weight = DataType_(1))
        {
            NodeType * node(new NodeType(id, weight));
            add_node(node);
        }
                
        /// sets the edgeweight bewteen nodes with the given nodeIDs to weight. If weight = 0, the edge  actually will be removed
        inline void set_edge_weight(int source_id, int destination_id, DataType_ weight)
        {
            set_edge_weight_internal(_node_mapping[source_id], _node_mapping[destination_id], weight);
        }

        /// sets the edgeweight between two nodes to weight
        inline void set_edge_weight(NodeType & source, NodeType & destination, DataType_ weight)
        {
            set_edge_weight(source.id(), destination.id(), weight);
        }

        /// adds an edge between nodes with the given nodeIDs and weight. In fact, if weight is 0, the edge will be removed
        inline void add_edge(int source_id, int destination_id, DataType_ weight = 1)
        {
            set_edge_weight(source_id, destination_id, weight);
        }

        /// adds an edge beween the nodes. 
        inline void add_edge(NodeType & source, NodeType & destination, DataType_ weight = 1)
        {
            set_edge_weight(source, destination, weight);
        }

        /// removes an probably existing edge between the nodes given by their IDs.
        inline void remove_edge(int source_id, int destination_id)
        {
            set_edge_weight(source_id, destination_id, 0);
        }

        /// removes an edge between the nodes
        inline void remove_edge(NodeType & source, NodeType & destination)
        {
            set_edge_weight(source, destination, 0);
        }

        /// returns the node at position index 
        inline NodeType & get_node(int index)
        {
            return *_nodes[index];
        }

        /// returns the node with the given ID, if it is part of this graph.
        inline NodeType &  get_node_by_id(int id)
        {
            return get_node(_node_mapping[id]);
        }

        /// returns the position index of a Node given by its id if it is part of this graph - or -1 otherwise.
        int get_node_index(int id)
        {
             if (_node_mapping.find(id) != _node_mapping.end())
                return _node_mapping[id];
            return -1;
        }

        /// the number of actually contained nodes.
        inline int node_count()
        {
            return _node_count;
        }
    };
}
#endif

