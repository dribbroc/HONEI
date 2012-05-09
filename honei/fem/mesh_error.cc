/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <markus.geveler@math.tu-dortmund.de>
 *
 * This file is part of the HONEI-FEM C++ library. LibLa is free software;
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

#include <honei/fem/mesh_error.hh>
#include <honei/util/stringify.hh>

#include <string>

using namespace honei;

MeshError::MeshError(const std::string & message) throw () :
    Exception(message)
{
}

MeshInternalIndexOutOfBounds::MeshInternalIndexOutOfBounds(unsigned long index, unsigned long max_index) throw () :
    MeshError("Internal index '" + stringify(index) + "' exceeds max index '"
            + stringify(max_index) + "'")
{
}
