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
        template<typename Tag_, typename MovementDirection_>
            class BoundaryInitFSI
            {
            };

        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_1> ///For solid approx. moving to direction 1, which means, the inter-/extrapolation direction is 3!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(Grid<D2Q9, DT_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_3)[packed_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_3)[pre_prev_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_3)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(Grid<D2Q9, DT_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_3)[packed_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8/15) * grid.d_x * solids.current_u) + (DT_(2/3) * v_m_1) - (DT_(2/5) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(Grid<D2Q9, DT_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_3)[packed_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] > grid.h->rows() * grid.h->columns() + 1) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8/15) * grid.d_y * solids.current_v) + (DT_(2/3) * v_m_1) - (DT_(2/5) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(Grid<D2Q9, DT_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_to_fluid_flags->lock(lm_read_only);

                            info.cuda_dir_3->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                (*data.h)[i] = (*solids.solid_to_fluid_flags)[i] ?
                                                    _extrapolation(grid, info, data, solids, i) :(*data.h)[i];

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !(*solids.solid_to_fluid_flags)[i]) ?
                                                    solids.current_u * grid.d_x : ((*solids.solid_to_fluid_flags)[i] ?
                                                            _interpolation_u(grid, info, data, solids, i) :(*data.u)[i]);

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !(*solids.solid_to_fluid_flags)[i]) ?
                                                    solids.current_v * grid.d_y : ((*solids.solid_to_fluid_flags)[i] ?
                                                            _interpolation_v(grid, info, data, solids, i) :(*data.v)[i]);
                            }

                            info.cuda_dir_3->unlock(lm_read_only);

                            solids.solid_to_fluid_flags->unlock(lm_read_only);
                            solids.boundary_flags->unlock(lm_read_only);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);
                        }
            };
    }
}
#endif
