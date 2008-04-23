/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBGRAPH_GUARD_GRAPH_ERROR_HH
#define LIBGRAPH_GUARD_GRAPH_ERROR_HH 1

#include <honei/util/exception.hh>

namespace honei
{
    /**
     * GraphError is the generic base class for all Graph related exceptions.
     *
     * \ingroup libgraph
     */
    class GraphError :
        public Exception
    {
        public:
            GraphError(const std::string & message) throw ();
    };

}

#endif

