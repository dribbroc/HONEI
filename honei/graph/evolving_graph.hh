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
#ifndef LIBGRAPH_GUARD_EVOLVING_GRAPH
#define LIBGRAPH_GUARD_EVOLVING_GRAPH 1

#include <tr1/memory>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/scale.hh>
#include <honei/la/sum.hh>
#include <honei/graph/abstract_graph.hh>
#include <honei/graph/graph.hh>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>
namespace honei
{

    template <typename DataType_> class EvolvingGraph:
        public AbstractGraph<DataType_>
    {
        private:
            typedef Node<DataType_> NodeType;
            typedef Graph<DataType_> GraphType;
            typedef DenseMatrix<DataType_> DM;
            typedef DenseVector<DataType_> DV;
            typedef SparseMatrix<DataType_> SM;
            std::vector<GraphType *> _slices;
            std::vector<DM *> _interpolation_coordinates;
//            std::map<int, NodeType> _nodes;
            std::vector<int> _slice_offset;
            
            DM _last_interpolation;
            float _last_alpha;;
           
            int _total_node_count;
            int _coordinate_dimensions;
            DataType_ _intertimeslice_weight;

            void assemble_coordinates()
            {
                this->_coordinates = new DM(_total_node_count, _coordinate_dimensions );
                for (int t(0); t < slice_count(); ++t)
                {
                
                    int offset = _slice_offset[t];
                    GraphType * slice = _slices[t];
                   
                    for (typename DM::ConstElementIterator i(slice->coordinates()->begin_elements()),
                            i_end(slice->coordinates()->end_elements()); i != i_end; ++i)
                    {
                        (*this->_coordinates)(offset + i.row(), i.column()) = *i;
                    }
                     
                }
            }

            /// creates and returns the complete nodeweight vector for the whole evolving graph
            void assemble_node_weights()
            {
                this->_node_weights = new DV(_total_node_count);
                for (int t(0); t < slice_count(); ++t)
                {
                    int offset = _slice_offset[t];
                    GraphType * slice = _slices[t];
                    for (typename DV::ConstElementIterator i(slice->node_weights()->begin_elements()),
                            i_end(slice->node_weights()->end_elements()); i != i_end; ++i)
                    {
                        (*this->_node_weights)[offset + i.index()] = *i;
                    }
                }
            }

            /// creates and returns the complete edge matrix for the whole evolviong graph.
void assemble_edges()
{
    this->_edges = new SM(_total_node_count, _total_node_count);
    for (int t(0); t < slice_count(); ++t)
    {
        int offset = _slice_offset[t];
        GraphType * slice = _slices[t];
        for (typename SM::NonZeroConstElementIterator i(slice->edges()->begin_non_zero_elements()),
                i_end(slice->edges()->end_non_zero_elements()); i != i_end; ++i)
        {
            (*this->_edges)(offset + i.row(), offset + i.column()) = *i;
        }
        if (t > 0)
        {
            GraphType * last_slice = _slices[t-1];
            for (int i = 0; i < slice->node_count(); ++i)
            {
                int last_index(last_slice->get_node_index(slice->get_node(i).id()));
                if (last_index >= 0)
                {
                    (*this->_edges)(_slice_offset[t-1] + last_index, offset + i) = _intertimeslice_weight;
                    (*this->_edges)(offset + i, _slice_offset[t-1] + last_index) = _intertimeslice_weight;
                }
            }
        }
    }
}


        public:
        /// creates an evolving graph with given number of dimensions for the node coordinates sets optionally the weight for intertimeslice edges
        EvolvingGraph(int coordinate_dimensions, DataType_ intertimeslice_weight = (DataType_)1) :
                _slices(),
//                _nodes(),
                _interpolation_coordinates(),
                _slice_offset(),
                _last_interpolation(1, 1),
                _total_node_count(0),
                _coordinate_dimensions(coordinate_dimensions),
                _intertimeslice_weight(intertimeslice_weight)
	{
	}

        virtual ~EvolvingGraph()
        {
          //  for(int i(0);  i < _interpolation_coordinates.size(); ++i)
          //      if (_interpolation_coordinates[i] != 0)
          //          delete(_interpolation_coordinates[i]);
        }

        /// adds a node to the evolving graph. this is necessary to put the same nodes to different timeslice-graphs
/*        void add_node(NodeType & node)
        {
            _nodes[node.id()] = node;
        } 
        
        /// adds a node to the evolving graph. this is necessary to put the same nodes to different timeslice-graphs
        void add_node(int timeslice_index, int node_id)
        {
            
            NodeType node(node_id);
            if (_nodes.find(node_id) == _nodes.end())
            {
                _nodes[node_id] = node;
            }
            else
                node = _nodes[node_id];            
            _slices[timeslice_index]->add_node(&node);
        } 
        
        void add_node(int timeslice_index, int node_id, DataType_ weight)
        {
            NodeType node(node_id, weight);
            if (_nodes.find(node_id) == _nodes.end())
            {
                _nodes[node_id] = node;                
            }
            else
            {
                node = _nodes[node_id];
                node.set_weight(weight);
            }        
            _slices[timeslice_index]->add_node(&node);
        } 
        

        /// returns the node with the given ID, if it was added.
        inline NodeType & get_node(int id)
        {
            return _nodes[id];
        }

        inline NodeType & get_node_by_id(int id)
        {
            return get_node(id);
        }
*/
        /// the number of DIFFERENT nodes contained in this graph - one node may occur many times in the assembled matrices 
        inline int node_count()
        {
            return _total_node_count;
        }
        

        /// the number of timeslices in this evolving graph
        inline int slice_count()
        {
            return _slices.size();
        }

        /// adds a timeslice to the evolving graph
        GraphType & add_timeslice(int node_count)
        {
            GraphType * g(new GraphType(node_count, 2));
            _slices.push_back(g);
            _slice_offset.push_back(_total_node_count);
            _total_node_count += node_count;
            std::cout << "added timeslice. total nodes: "<< _total_node_count << " (" << node_count << "new)\n";
            return *_slices.back();
        }

        /// creates and returns the complete coordinate matrix for the whole evolving graph

        /// returns the index of the timeslice that contains the given nodeIndex
        int timeslice_index(int node_index)
        {
            int t(0);
            for (; t < slice_count()  - 1; ++t)
            {
                if (node_index < _slice_offset[t+1])
                    return t;
            }
            return slice_count()-1;
        }

        /// returns true if the indices are in the same timeslice - false otherwise
        virtual bool same_timeslice(int index1, int index2)
        {
            int slice_index = timeslice_index(index1);
            return (slice_index == slice_count()-1 || index2 < _slice_offset[slice_index+1]) &&
                index2 >= _slice_offset[slice_index];
        }

        /// returns the timslice graph at a given time t
        inline GraphType & get_timeslice(int t)
        {
            return *_slices[t];
        }
        
        inline GraphType & last_timeslice()
        {
            return *_slices.back();
        }

        DenseMatrix<DataType_> * coordinates()
        {
            if (!this->_coordinates)
                assemble_coordinates();
            return this->_coordinates;
        }

        DenseVector<DataType_> * node_weights()
        {
            if (!this->_node_weights)
                assemble_node_weights();
            return this->_node_weights;
        }

        SparseMatrix<DataType_> * edges()
        {
            if (!this->_edges)
                assemble_edges();
            return this->_edges;
        }

        void reassemble_graph()
        {
            std::cout <<"reassemble graph\n";
           if(this->_coordinates)
              delete(this->_coordinates);
            std::cout <<"reassemble coords\n";
            assemble_coordinates();
            
            if(this->_node_weights)
               delete(this->_node_weights);
            std::cout <<"reassemble weights\n";
            assemble_node_weights();
            if(this->_edges)
                delete(this->_edges);
            assemble_edges();
        }
        
        void update_slice_coordinates(DM & new_coordinates)
        {
       //     std::cout << "update slice coordinates\n";
            for (int timeslice_idx(0); timeslice_idx < slice_count(); ++timeslice_idx)
            {
                DM * slice_coordinates(_slices[timeslice_idx]->coordinates());
                typename DM::ConstElementIterator k(new_coordinates.begin_elements());
                k += _slice_offset[timeslice_idx] * _coordinate_dimensions;
            
                for (typename DM::ElementIterator i(slice_coordinates->begin_elements()), end_i(slice_coordinates->end_elements());
                i != end_i; ++i, ++k)
                {                
                    *i = *k;
                }
            }                
        }
        
        void update_slice_coordinates()
        {
      //      std::cout << "update slice coordinates\n";
            for (int timeslice_idx(0); timeslice_idx < slice_count(); ++timeslice_idx)
            {
                DM * slice_coordinates(_slices[timeslice_idx]->coordinates());
                typename DM::ConstElementIterator k(this->_coordinates->begin_elements());
                k += _slice_offset[timeslice_idx] * _coordinate_dimensions;
            
                for (typename DM::ElementIterator i(slice_coordinates->begin_elements()), end_i(slice_coordinates->end_elements());
                i != end_i; ++i, ++k)
                {                
                    *i = *k;
                }
            }                
        }
        
        virtual bool includes_timeslices()
        {
            return true;
        }
    };
}
#endif
