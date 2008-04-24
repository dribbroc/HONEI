/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LibLBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBLBM_GUARD_TAGS_HH
#define LIBLBM_GUARD_TAGS_HH 1

/**
 * \file
 * Implementation of a SWE solver using LBM.
 *
 * \ingroup grpliblbm
 **/

namespace honei
{

    namespace lbm
    {
        namespace lbm_lattice_types
        {
            class D2Q9;
        }
        namespace lbm_grid_types
        {
            class RECTANGULAR;
        }
        namespace lbm_source_types
        {
            class SIMPLE;
        }
        namespace lbm_source_schemes
        {
            class BASIC;
            class CENTERED;
        }
        namespace lbm_boundary_types
        {
            class NOSLIP_PERIODIC;
        }

        namespace lbm_applications
        {
            class LABSWE;
        }

        namespace lbm_directions
        {
            class X;
            class Y;
        }
    }

}
#endif
