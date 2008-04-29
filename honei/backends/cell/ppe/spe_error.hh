/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'exception.hh' from Paludis, which is:
 *     Copyright (c) 2005, 2006, 2007 Ciaran McCreesh
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef BACKENDS_GUARD_CELL_PPE_SPE_ERROR_HH
#define BACKENDS_GUARD_CELL_PPE_SPE_ERROR_HH 1

#include <honei/util/exception.hh>

#include <string>

namespace honei
{
    /**
     * SPEError is thrown by SPEManager and related classes whenever an error
     * occurs in interfacing Libspe2.
     *
     * \ingroup grpexceptions
     * \ingroup grpcell
     */
    struct SPEError :
        public ExternalError
    {
        /**
         * Constructor.
         *
         * \param msg The error message.
         * \param reason The reason for the error message.
         */
        SPEError(const std::string & function, int errno) throw ();
    };
}

#endif
