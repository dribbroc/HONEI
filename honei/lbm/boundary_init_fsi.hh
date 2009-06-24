/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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

#ifndef LBM_GUARD_BOUNDARY_INIT_FSI_HH
#define LBM_GUARD_BOUNDARY_INIT_FSI_HH 1

namespace honei
{
    namespace lbm
    {
        template<typename Tag_>
            class BoundaryInitFSI
            {
            };

        template<>
            class BoundaryInitFSI<tags::CPU>
            {
                public:
                    template<typename DT_>
                        static void value(Grid<D2Q9, DT1_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT1_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_to_fluid_flags->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !(*solids.solid_to_fluid_flags)[i]) ?
                                                    solids.current_u * grid.d_x :
                                                        (*data.u)[i];

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !(*solids.solid_to_fluid_flags)[i]) ?
                                                    solids.current_v * grid.d_y :
                                                        (*data.v)[i];
                            }

                            solids.solid_to_fluid_flags->unlock(lm_read_only);
                            solids.boundary_flags->unlock(lm_read_only);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);
                        }
            };
    }
}
#endif
