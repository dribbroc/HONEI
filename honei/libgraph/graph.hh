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

#ifndef LIBGRAPH_GUARD_GRAPH
#define LIBGRAPH_GUARD_GRAPH 1

#include <tr1/memory>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libgraph/node.hh>
#include <honei/libgraph/abstract_graph.hh>
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
        int _nodeCount;
        int _maxNodes;
        int _coordinateDimensions;
    std::map<int, int> _nodeMapping;

        /// sets edges[v_index][w_index] = edges[w_index][v_index] = weight
        void setEdgeWeightInternal(int v_index, int w_index, DataType_ weight)
        {
            (*this->_edges)(v_index, w_index) = weight;
            (*this->_edges)(w_index, v_index) = weight;
        }

    public:
        /// constructs new graph with given number of nodes and a given dimension for the coordinates
        Graph(int nodes, int coordinateDimensions) :
            _coordinateDimensions(coordinateDimensions),
            _nodeMapping(),
            _nodeCount(0),
            _maxNodes(nodes)
        {
            _nodes = new NodeType*[nodes];
        this->_coordinates = new DenseMatrix<DataType_>(nodes, coordinateDimensions);
            this->_edges = new SparseMatrix<DataType_>(nodes, nodes);
            this->_nodeWeights = new DenseVector<DataType_>(nodes, (DataType_)1);
        }

        /// adds a node to this graph, puts its initial position and its weight into the relevant matrices
        void addNode(NodeType * node)
        {
            if (_nodeCount < _maxNodes)
            {
                _nodes[_nodeCount] = node;
                _nodeMapping[node->getID()] = _nodeCount;
                for (int i(0); i < _coordinateDimensions; ++i)
                {
                    (*this->_coordinates)(_nodeCount, i) = 2 - (DataType_)std::rand() / RAND_MAX * 4; // i < node->getPosition()->size() ? (*node->getPosition())[i] : std::rand();
                }
                (*this->_nodeWeights)[_nodeCount] = node->getWeight();
                ++_nodeCount;
            }
        }

        /// sets the edgeweight bewteen nodes with the given nodeIDs to weight. If weight = 0, the edge  actually will be removed
        inline void setEdgeWeight(int sourceID, int destinationID, DataType_ weight)
        {
            setEdgeWeightInternal(_nodeMapping[sourceID], _nodeMapping[destinationID], weight);
        }

        /// sets the edgeweight between two nodes to weight
        inline void setEdgeWeight(NodeType * source, NodeType * destination, DataType_ weight)
        {
            setEdgeWeight(source->getID(), destination->getID(), weight);
        }

        /// adds an edge between nodes with the given nodeIDs and weight. In fact, if weight is 0, the edge will be removed
        inline void addEdge(int sourceID, int destinationID, DataType_ weight = 1)
        {
            setEdgeWeight(sourceID, destinationID, weight);
        }

        /// adds an edge beween the nodes. 
        inline void addEdge(NodeType * source, NodeType * destination, DataType_ weight = 1)
        {
            setEdgeWeight(source, destination, weight);
        }

        /// removes an probably existing edge between the nodes given by their IDs.
        inline void removeEdge(int sourceID, int destinationID)
        {
            setEdgeWeight(sourceID, destinationID, 0);
        }

        /// removes an edge between the nodes
        inline void removeEdge(NodeType * source, NodeType * destination)
        {
            setEdgeWeight(source, destination, 0);
        }

        /// returns the node at position index 
        inline NodeType * getNode(int index)
        {
            return _nodes[index];
        }

        /// returns the node with the given ID, if it is part of this graph.
        inline NodeType *  getNodeByID(int ID)
        {
            return getNode(_nodeMapping[ID]);
        }

        /// returns the position index of a Node given by its id if it is part of this graph - or -1 otherwise.
        int getNodeIndex(int id)
        {
             if (_nodeMapping.find(id) != _nodeMapping.end())
                return _nodeMapping[id];
            return -1;
        }

        /// the number of actually contained nodes.
        inline int nodeCount()
        {
            return _nodeCount;
        }
    };
}
#endif
