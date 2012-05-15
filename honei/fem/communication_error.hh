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

#pragma once
#ifndef FEM_GUARD_COMMUNICATION_ERROR_HH
#define FEM_GUARD_COMMUNICATION_ERROR_HH 1

#include <honei/util/exception.hh>

#include <string>

namespace honei
{
    class CommunicationError :
        public Exception
    {
        public:
            CommunicationError(const std::string & message) throw ();
    };

    class CommunicationHaloOverlapMismatch :
        public CommunicationError
    {
        public:
            CommunicationHaloOverlapMismatch(unsigned long comm_overlap, unsigned long halo_overlap) throw ();
    };
}

#endif
