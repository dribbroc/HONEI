/* vim: set sw=4 sts=4 et nofoldenable : */

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

#ifndef LIBGRAPH_GUARD_NODE
#define LIBGRAPH_GUARD_NODE 1

#include <iostream>
#include <tr1/memory>
#include <honei/libla/dense_vector.hh>

namespace honei
{
    template <typename DataType_> class Node
    {
    private:
        // Graph<DataType_> _graph*;        // the graph this node belongs to
        unsigned long _id;            // key to represent the node
        DataType_ _weight;    // initial weight of the node
        typedef DenseVector<DataType_> DV; // initial position of this node
        DV *  _position;

    public:
        /// creates a new node with a given id, weight and a position defined by a DenseVector
        Node(unsigned long ID, DataType_ weight, DV * position) :
                    _position(position),
                    _id(ID),
                    _weight(weight)
        {
        }

        Node (unsigned long ID, DataType_ weight):
            _position(new DV(2)),
            _id(ID),
            _weight(weight)
        {
        }

        /// creates a new node with a given id, weight and coordinates (x,y).
        Node(unsigned long ID, DataType_ weight, DataType_ x, DataType_ y) :
            _position(new DV(2)),
            _id(ID),
            _weight(weight)
        {
            (*_position)[0] = x;
            (*_position)[1] = y;
        }

        /// creates a new node with given id, weight and coordnates (x,y,z) - for people who like 3D graphs^^
        Node(unsigned long ID, DataType_ weight, DataType_ x, DataType_ y, DataType_ z) :
            _position(new DV(3)),
            _id(ID),
            _weight(weight)
        {
            (*_position)[0] = x;
            (*_position)[1] = y;
            (*_position)[2] = z;
        }

        ~Node()
        {
            if (_position)
                delete(_position);
        }

        /// returns the initial weight of this node
        inline DataType_ getWeight()
        {
            return _weight;
        }

        /// sets the initial weight, which is put into graph's nodeWeight vector while adding the node
        inline void setWeight(DataType_ weight)
        {
                _weight = weight;
        }

        /// returns the all important node ID
        inline unsigned long getID()
        {
            return _id;
        }

        /// returns the initial position of this node.
        inline DenseVector<DataType_> * getPosition()
        {
            return _position;
        }

        /// sets the initial position.it i s put into graph's coordinate matrix when adding.
        inline void setPosition(DenseVector<DataType_> * position)
        {
            _position = position;
        }
    };
}
#endif
