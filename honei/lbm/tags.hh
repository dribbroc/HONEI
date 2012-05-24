/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#pragma once
#ifndef LBM_GUARD_TAGS_HH
#define LBM_GUARD_TAGS_HH 1

/**
 * \file
 * Definition of LBM related tags.
 *
 * \ingroup grpliblbm
 **/

namespace honei
{

    namespace lbm
    {
        namespace lbm_lattice_types
        {
            class D2Q9
            {
                public:
                    class DIR_ODD;
                    class DIR_EVEN;

                    class DIR_0;
                    class DIR_1;
                    class DIR_2;
                    class DIR_3;
                    class DIR_4;
                    class DIR_5;
                    class DIR_6;
                    class DIR_7;
                    class DIR_8;
            };
        }
        namespace lbm_grid_types
        {
            class RECTANGULAR;
        }
        namespace lbm_force
        {
            class SIMPLE;
            class CONSTANT;
            class CENTRED;
            class NONE;
        }
        namespace lbm_source_schemes
        {
            class BASIC;
            class BED_SLOPE;
            class BED_FRICTION;
            class BED_FULL;
            class NONE;
        }
        namespace lbm_boundary_types
        {
            class NOSLIP_PERIODIC;
            class NOSLIP;
            class NOSLIP_2;
            class PERIODIC;
            class DRIVEN_CAVITY;
            class DIRICHLET_SLIP;
        }

        namespace lbm_applications
        {
            class LABSWE;
            class LABNAVSTO;
        }

        namespace lbm_directions
        {
            class X;
            class Y;
        }

        namespace lbm_modes
        {
            class DRY;
            class WET;
        }
    }

}
#endif
