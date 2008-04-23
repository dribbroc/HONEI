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
#include <honei/la/dense_vector.hh>

namespace honei
{
    template <typename DataType_> class Node
    {
    private:
        unsigned long _id;            // key to represent the node
        DataType_ _weight;    // initial weight of the node


    public:

        Node () :
            _id(0),
            _weight(1)
        {
        }

        Node (unsigned long ID, DataType_ weight = 1):
            _id(ID),
            _weight(weight)
        {
        }


        /// returns the initial weight of this node
        inline DataType_ get_weight()
        {
            return _weight;
        }

        /// sets the initial weight, which is put into graph's nodeWeight vector while adding the node
        inline void set_weight(DataType_ weight)
        {
                _weight = weight;
        }

        /// returns the all important node ID
        inline unsigned long id()
        {
            return _id;
        }
    };
}
#endif
