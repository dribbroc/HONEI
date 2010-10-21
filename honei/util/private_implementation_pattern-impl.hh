/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Base upon 'private_implementation_pattern-impl.hh' from Paludis, which is:
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
#ifndef LIBUTIL_GUARD_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH
#define LIBUTIL_GUARD_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH 1

#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    template <typename T_>
    PrivateImplementationPattern<T_, Shared>::PrivateImplementationPattern(Implementation<T_> * imp) :
        _imp(imp)
    {
    }

    template <typename T_>
    PrivateImplementationPattern<T_, Shared>::PrivateImplementationPattern(shared_ptr<Implementation<T_> > imp) :
        _imp(imp)
    {
    }

    template <typename T_>
    PrivateImplementationPattern<T_, Shared>::~PrivateImplementationPattern()
    {
    }

    template <typename T_>
    PrivateImplementationPattern<T_, Single>::PrivateImplementationPattern(Implementation<T_> * imp) :
        _imp(imp)
    {
    }

    template <typename T_>
    PrivateImplementationPattern<T_, Single>::~PrivateImplementationPattern()
    {
        delete _imp;
    }
}

#endif
