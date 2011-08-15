/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#pragma once
#ifndef LIBMATH_GUARD_METHODS_HH
#define LIBMATH_GUARD_METHODS_HH 1
/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

namespace methods
{
    struct CG
    {
    };

    struct JAC
    {
    };

    struct SPAI
    {
    };

    struct ILU
    {
    };

    struct NONE
    {
    };

    struct VAR
    {
    };

    struct PCG
    {
        public:
            struct JAC;
    };

    struct CYCLE
    {
        public:
            struct F;
            struct W
            {
                struct STATIC;
            };
            struct V
            {
                struct STATIC;
            };
    };
    struct FIXED;
    struct MIXED;

    struct PROLMAT;

    struct NATURAL;
    struct TWO_LEVEL;
}


namespace applications
{
    struct POISSON;
}
namespace boundary_types
{
    struct DIRICHLET
    {
        public:
            struct DIRICHLET_0;
    };
    struct NEUMANN;
    struct DIRICHLET_NEUMANN
    {
    };
}

#endif
