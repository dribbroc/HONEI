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
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libgraph/node.hh>
#include <map>

namespace honei 
{
    template <typename DataType_> 
    class AbstractGraph
    {
    protected:
        DenseMatrix<DataType_> * _coordinates;
        SparseMatrix< DataType_ > * _edges;
        DenseVector<DataType_> * _nodeWeights;
    public:
        virtual void addNode(Node<DataType_> * Node) = 0;

    virtual Node<DataType_> * getNode(int index) = 0;

    virtual Node<DataType_> * getNodeByID(int id) = 0;

    AbstractGraph() :        
        _coordinates(0),
        _edges(0),
        _nodeWeights(0)
    {        
    }

    virtual DenseMatrix<DataType_> * coordinates()
    {
        return _coordinates;
    }

    virtual DenseVector<DataType_> * nodeWeights()
    {
        return _nodeWeights;
    }
    
    virtual SparseMatrix<DataType_> * edges()
    {
        return _edges;
    }

    virtual int nodeCount() = 0;
        
    virtual bool sameTimeslice(int index1, int index2)
    {
        return true;
    }
    };
}

#endif
