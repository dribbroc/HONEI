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
#ifndef LIBGRAPH_GUARD_EVOLVING_ANIMATOR
#define LIBGRAPH_GUARD_EVOLVING_ANIMATOR 1

#include <tr1/memory>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/scale.hh>
#include <honei/la/sum.hh>
#include <honei/graph/evolving_graph.hh>
#include <honei/graph/graph.hh>
#include <map>
#include <iostream>
#include <math.h>
namespace honei
{

    template <typename Tag_, typename DataType_> 
    class EvolvingAnimator
    {
        private:
            typedef Node<DataType_> NodeType;
            typedef Graph<DataType_> GraphType;
            typedef DenseMatrix<DataType_> DM;
            typedef DenseVector<DataType_> DV;
            typedef SparseMatrix<DataType_> SM;
            
            EvolvingGraph<DataType_> & _graph;
            float _step_size;
            int _timeslice_idx;
            float _alpha;
            float _last_time;     
            
            std::vector<DM *> _interpolation_coordinates;            
            DM _last_interpolation;     

            

        public:
        /// creates an evolving graph with given number of dimensions for the node coordinates sets optionally the weight for intertimeslice edges
        EvolvingAnimator(EvolvingGraph<DataType_> & graph, float step_size):
            _graph(graph),
            _step_size(step_size),
            _timeslice_idx(0),
            _alpha(0.0f),
            _last_time(-1.0f),
            _last_interpolation(1, 1)
        {
            
        }
        
        ~EvolvingAnimator()
        {
          //  for(int i(0);  i < _interpolation_coordinates.size(); ++i)
          //      if (_interpolation_coordinates[i] != 0)
          //          delete(_interpolation_coordinates[i]);
        }

       

        /// returns the timslice graph at a given time t
        inline int timeslice_index()
        {
            return _timeslice_idx >= _graph.slice_count() ? _graph.slice_count() - 1 : _timeslice_idx;
        }
        
        inline float time()
        {
            return _last_time;
        }

        inline DenseMatrix<DataType_> & coordinates()
        {
            return _last_interpolation;
        }

        inline DenseVector<DataType_> * node_weights()
        {
            return _graph.get_timeslice(timeslice_index()).node_weights();
        }

        inline SparseMatrix<DataType_> * edges()
        {
            return _graph.get_timeslice(timeslice_index()).edges();
        }
        
        void prepare_interpolation()
        {            
            _last_time = 0.0f;
            _last_interpolation = _graph.get_timeslice(0).coordinates()->copy();
            
            for(unsigned long i(0);  i < _interpolation_coordinates.size(); ++i)
                if (_interpolation_coordinates[i] != 0)
                    delete(_interpolation_coordinates[i]);
            _interpolation_coordinates.clear();
            
            
            
            for (int slice_idx(0); slice_idx < _graph.slice_count()-1; ++slice_idx)
            {
                DM * slice_begin_coordinates(_graph.get_timeslice(slice_idx).coordinates());
                DM * next_slice_coordinates(_graph.get_timeslice(slice_idx+1).coordinates());
                DM * slice_final_coordinates = new DM(slice_begin_coordinates->rows(), slice_begin_coordinates->columns());
                       
                for (typename DM::ElementIterator i(slice_final_coordinates->begin_elements()), i_end(slice_final_coordinates->end_elements()),
                    k(slice_begin_coordinates->begin_elements()); i != i_end; ++i, ++k)
                {
                    NodeType & node(_graph.get_timeslice(slice_idx).get_node(i.row()));
                    int corresponding_idx_next_slice = _graph.get_timeslice(slice_idx+1).get_node_index(node.id());
                    *i = corresponding_idx_next_slice == -1 ? *k : (*next_slice_coordinates)[corresponding_idx_next_slice][i.column()];
                }                
                _interpolation_coordinates.push_back(slice_final_coordinates);
            }
                            
            _interpolation_coordinates.push_back(_graph.last_timeslice().coordinates());
        }
        
        inline DM next_step()
        {
            return interpolate_timeslice(_last_time + _step_size);        
        }
        
        inline DM rewind()
        {
            return interpolate_timeslice(0.0f);
        }
        
        inline bool end()
        {
            return _last_time >= _graph.slice_count();
        }
        
        //template <typename Tag_>
        DM interpolate_timeslice(float time)
        {
            if (time == _last_time)
                return _last_interpolation;
            _last_time = time;
            
            _timeslice_idx = (int)floor(time);
            _alpha = time - _timeslice_idx;
            if (end())
            {
                _last_interpolation = *_graph.last_timeslice().coordinates();
            /*    _timeslice_idx = _graph.slice_count() - 1;
                _alpha = 0;
                _last_time =  _timeslice_idx + _step_size;
            */
            }
            else
            {
                DM * coordinates_begin(_graph.get_timeslice(_timeslice_idx).coordinates());
                DM * coordinates_final(_interpolation_coordinates[_timeslice_idx]);
                
                _last_interpolation = coordinates_begin->copy();
                Scale<Tag_>::value(_last_interpolation, 1.0f - _alpha);
                DM tmp = coordinates_final->copy();
                Sum<Tag_>::value(_last_interpolation, Scale<Tag_>::value(tmp, _alpha));
            }
            return _last_interpolation;           
        }
    };
}
#endif
