/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Base upon 'private_implementation_pattern.hh' from Paludis, which is:
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

#pragma once
#ifndef LIBUTIL_GUARD_PRIVATE_IMPLEMENTATION_PATTERN_HH
#define LIBUTIL_GUARD_PRIVATE_IMPLEMENTATION_PATTERN_HH 1

#include <memory>

namespace honei
{
    /**
     * \{
     * \name Private implementation pattern storage tags
     *
     * Possible tags that govern the storage of private implementations.
     */

    /**
     * Tag for a shared private implementation.
     */
    struct Shared;

    /**
     * Tag for a single private implementation.
     */
    struct Single;

    /// \}

    /**
     * Private implementation data, to be specialised for any class that uses
     * PrivateImplementationPattern.
     */
    template <typename T_> struct Implementation;

    /**
     * \{
     *
     * PrivateImplementationPattern is a utility base class that can be used to
     * easily create classes with private implementation.
     */

    template <typename T_, typename Method_> class PrivateImplementationPattern;

    template <typename T_> class PrivateImplementationPattern<T_, Shared>
    {
        protected:
            /// Our implementation.
            std::shared_ptr<Implementation<T_> > _imp;

        public:
            explicit PrivateImplementationPattern(Implementation<T_> * imp);

            explicit PrivateImplementationPattern(std::shared_ptr<Implementation<T_> > imp);

            ~PrivateImplementationPattern();
    };

    template <typename T_> class PrivateImplementationPattern<T_, Single>
    {
        protected:
            /// Our implementation.
            Implementation<T_> * _imp;

        public:
            explicit PrivateImplementationPattern(Implementation<T_> * imp);

            ~PrivateImplementationPattern();
    };

    /// \}
}

#endif
