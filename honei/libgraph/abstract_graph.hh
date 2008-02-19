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

#ifndef LIBGRAPH_GUARD_ABSTRACT_GRAPH
#define LIBGRAPH_GUARD_ABSTRACT_GRAPH 1

#include <tr1/memory>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libgraph/node.hh>
#include <map>

namespace honei 
{
    template <typename DataType_> 
    class AbstractGraph
    {
    protected:
        DenseMatrix<DataType_> * _coordinates;
        SparseMatrix< DataType_ > * _edges;
        DenseVector<DataType_> * _node_weights;
    public:
        AbstractGraph() :
            _coordinates(0),
            _edges(0),
            _node_weights(0)
        {
        }

        ~AbstractGraph()
        {
            delete _coordinates;
            delete _edges;
            delete _node_weights;
        }

        virtual inline DenseMatrix<DataType_> * coordinates()
        {
            return _coordinates;
        }

        virtual inline DenseVector<DataType_> * node_weights()
        {
            return _node_weights;
        }

        virtual inline SparseMatrix<DataType_> * edges()
        {
            return _edges;
        }

        virtual int node_count() = 0;
                
        virtual inline int timeslice_index(int node_index)
        {
            return 0;
        }

        virtual inline bool same_timeslice(int index1, int index2)
        {
            return true;
        }
    };
}

#endif

