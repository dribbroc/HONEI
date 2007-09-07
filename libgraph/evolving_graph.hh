/* vim: set sw=4 sts=4 et nofoldenable syntax on */

/*
 * Copyright (c) 2007 Thorsten Deinert <thosten.deinert@uni-dortmund.de>
 *
 * This file is part of the Graph C library. LibGraph is free software;
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

#include <iostream>
#include <tr1/memory>
#include <iostream>
#include <cmath>
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libgraph/graph.hh>
#include <map>

namespace honei
{

    template <typename DataType_> class EvolvingGraph
    {
    private:
        typedef Node<DataType_> NodeType;
        typedef Graph<DataType_> GraphType;
        typedef DenseMatrix<DataType_> DM;
        typedef DenseVector<DataType_> DV;
        typedef SparseMatrix<DataType_> SM;
        std::vector<GraphType *> _slices;
        std::map<int, NodeType *> _nodes;
        std::vector<int> _sliceOffset;
        int _totalNodeCount;
        int _coordinateDimensions;
        DataType_ _interTimesliceWeight;
    public:
        /// creates an evolving graph with given number of dimensions for the node coordinates sets optionally the weight for intertimeslice edges
        EvolvingGraph(int coordinateDimensions, DataType_ interTimesliceWeight = 1) :
            _slices(),
            _nodes(),
            _sliceOffset(),
            _totalNodeCount(0),
            _coordinateDimensions(coordinateDimensions),
            _interTimesliceWeight(interTimesliceWeight)
        {
        }

        /// adds a node to the evolving graph. this is necessary to put the same nodes to different timeslice-graphs
        void addNode(NodeType * node)
        {
            _nodes[node->getID()] = node;
        }

        /// returns the node with the given ID, if it was added.
        NodeType * getNode(int id)
        {
            return _nodes[id];
        }

        /// the number of nodes contained
        int nodeCount()
        {
            return _nodes.size();
        }

        /// the number of timeslices in this evolving graph
        int sliceCount()
        {
            return _slices.size();
        }

        /// adds a timeslice to the evolving graph
        void addTimeslice(GraphType * sliceGraph)
        {
            _slices.push_back(sliceGraph);
            _sliceOffset.push_back(_totalNodeCount);
            _totalNodeCount += sliceGraph->nodeCount();
        }

        /// creates and returns the complete coordinate matrix for the whole evolving graph
        DM * getCoordinates()
        {
            std::cout << "Total Nodes: " << _totalNodeCount << "\n";
            DM * coordinates = new DM(_totalNodeCount, _coordinateDimensions);
            for (int t(0); t < sliceCount(); ++t)
            {
                int offset = _sliceOffset[t];
                GraphType * slice = _slices[t];
                for (typename DM::ConstElementIterator i(slice->getCoordinates()->begin_elements()),
                    i_end(slice->getCoordinates()->end_elements()); i != i_end; ++i)
                {
                    (*coordinates)[i.row()][offset + i.column()] = *i;
                }
            }
            return coordinates;
        }

        /// creates and returns the complete nodeweight vector for the whole evolving graph
        DV * getNodeWeights()
        {
            DV * nodeWeights = new DV(_totalNodeCount);

            for (int t(0); t < sliceCount(); ++t)
            {
                int offset = _sliceOffset[t];
                GraphType * slice = _slices[t];
                for (typename DV::ConstElementIterator i(slice->getNodeWeights()->begin_elements()),
                    i_end(slice->getNodeWeights()->end_elements()); i != i_end; ++i)
                {
                    (*nodeWeights)[offset + i.index()] = *i;
                }
            }

            return nodeWeights;
        }

        /// creates and returns the complete edge matrix for the whole evolviong graph.
        SM * getEdges()
        {
            SM * edgeWeights = new SM(_totalNodeCount, _totalNodeCount);
            for (int t(0); t < sliceCount(); ++t)
            {
                int offset = _sliceOffset[t];
                GraphType * slice = _slices[t];
                for (typename SM::ConstElementIterator i(slice->getEdges()->begin_non_zero_elements()),
                    i_end(slice->getEdges()->end_non_zero_elements()); i != i_end; ++i)
                {
                    (*edgeWeights)[offset + i.row()][offset + i.column()] = *i;
                }
                if (t > 0)
                {
                    GraphType * last_slice = _slices[t-1];
                    for (int i = 0; i < slice->nodeCount(); ++i)
                    {
                        int last_index(last_slice->getNodeIndex(slice->getNode(i)->getID()));
                        if (last_index >= 0)
                        {
                            (*edgeWeights)[_sliceOffset[t-1] + last_index][offset + i] = _interTimesliceWeight;
                            (*edgeWeights)[offset + i][_sliceOffset[t-1] + last_index] = _interTimesliceWeight;
                        }
                    }
                }
            }
            return edgeWeights;
        }

        /// returns the index of the timeslice that contains the given nodeIndex
        int getTimesliceIndex(int nodeIndex)
        {
            int t(0);
            for (; t < sliceCount()  - 1; ++t)
            {
                if (nodeIndex < _sliceOffset[t+1])
                return t;
            }
            return t+1;
        }

        /// returns true if the indices are in the same timeslice - false otherwise
        bool SameTimeslice(int index1, int index2)
        {
            int sliceIndex = getTimesliceIndex(index1);
            return     (sliceIndex == sliceCount()-1 || index2 < _sliceOffset[sliceIndex+1]) &&
                    index2 >= _sliceOffset[sliceIndex];
        }

        /// returns the timslice graph at a given time t 
        GraphType * getTimeslice(int t)
        {
            return _slices[t];
        }
    };
}