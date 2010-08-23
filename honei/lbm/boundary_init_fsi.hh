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

#pragma once
#ifndef LBM_GUARD_BOUNDARY_INIT_FSI_HH
#define LBM_GUARD_BOUNDARY_INIT_FSI_HH 1

#include <climits>
#include <honei/lbm/tags.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/grid.hh>

namespace honei
{
    namespace lbm
    {
        template<typename Tag_, typename MovementDirection_>
            class BoundaryInitFSI
            {
            };

        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_1> ///For solid approx. moving to direction 1, which means, the inter-/extrapolation direction is 5!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_5)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_5)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_5)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_5)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_5)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_5)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_5)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_5)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_5)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_5)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_5)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_5)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_5)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_5)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_5)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_5)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_5)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_5)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_5)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_5)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_5)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_5->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_5->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };

        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_1>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_2> ///For solid approx. moving to direction 2, which means, the inter-/extrapolation direction is 6!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_6)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_6)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_6)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_6)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_6)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_6)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_6)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_6)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_6)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_6)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_6)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_6)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_6)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_6)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_6)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_6)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_6)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_6)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_6)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_6)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_6)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_6->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_6->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_3> ///For solid approx. moving to direction 2, which means, the inter-/extrapolation direction is 7!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_7)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_7)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_7)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_7)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_7)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_7)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_7)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_7)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_7)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_7)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_7)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_7)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_7)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_7)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_7)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_7)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_7)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_7)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_7)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_7)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_7)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_7->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_7->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_4> ///For solid approx. moving to direction 4, which means, the inter-/extrapolation direction is 8!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_8)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_8)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_8)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_8)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_8)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_8)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_8)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_8)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_8)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_8)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_8)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_8)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_8)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_8)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_8)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_8)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_8)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_8)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_8)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_8)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_8)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_8->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_8->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_5> ///For solid approx. moving to direction 5, which means, the inter-/extrapolation direction is 1!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_1)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_1)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_1)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_1)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_1)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_1)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_1)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_1)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_1)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_1)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_1)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_1)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_1)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_1)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_1)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_1)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_1)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_1)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_1)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_1)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_1)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_1->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_1->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_6> ///For solid approx. moving to direction 6, which means, the inter-/extrapolation direction is 2!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_2)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_2)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_2)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_2)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_2)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_2)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_2)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_2)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_2)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_2)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_2)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_2)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_2)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_2)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_2)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_2)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_2)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_2)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_2)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_2)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_2)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_2->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_2->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_7> ///For solid approx. moving to direction 7, which means, the inter-/extrapolation direction is 3!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_3)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_3)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_3)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_3)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_3)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_3)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_3)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_3)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_3)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_3)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_3->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_3->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template<>
            class BoundaryInitFSI<tags::CPU, D2Q9::DIR_8> ///For solid approx. moving to direction 8, which means, the inter-/extrapolation direction is 4!
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {
                            ///Determine extrapolation_indices:
                            bool prev(!((*info.cuda_dir_4)[packed_index] ==  ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_4)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_4)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_4)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_4)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_4)[prev_index] : prev_index);

                            bool pre_pre_prev(pre_prev ? (!((*info.cuda_dir_4)[pre_prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_4)[pre_prev_index]] ? true : false) : false) : false);
                            unsigned long pre_pre_prev_index(pre_pre_prev ? (*info.cuda_dir_4)[pre_prev_index] : pre_prev_index);

                            ///Gather extrapolation values:
                            DT_ v_m_1((*data.h)[prev_index]);
                            DT_ v_m_2((*data.h)[pre_prev_index]);
                            DT_ v_m_3((*data.h)[pre_pre_prev_index]);

                            return (DT_(3.) * (v_m_1 - v_m_2) + v_m_3);

                        }

                    template<typename DT_>
                        static DT_ _interpolation_u(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_4)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_4)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_4)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_4)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_4)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_4)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.u)[prev_index]);
                            DT_ v_m_2((*data.u)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_u) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }

                    template<typename DT_>
                        static DT_ _interpolation_v(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          unsigned long packed_index)
                        {

                            bool prev(!((*info.cuda_dir_4)[packed_index] == ULONG_MAX) ?
                                    (!(*solids.solid_flags)[(*info.cuda_dir_4)[packed_index]] ? true : false) : false);
                            unsigned long prev_index(prev ? (*info.cuda_dir_4)[packed_index] : packed_index);

                            bool pre_prev(prev ? (!((*info.cuda_dir_4)[prev_index] == ULONG_MAX) ?
                                        (!(*solids.solid_flags)[(*info.cuda_dir_4)[prev_index]] ? true : false) : false) : false);
                            unsigned long pre_prev_index(pre_prev ? (*info.cuda_dir_4)[prev_index] : prev_index);

                            //assuming q = 1/2
                            DT_ v_m_1((*data.v)[prev_index]);
                            DT_ v_m_2((*data.v)[pre_prev_index]);

                            return (DT_(8./15.) * solids.current_v) + (DT_(2./3.) * v_m_1) - (DT_(2./5.) * v_m_2);
                        }
                public:
                    template<typename DT_>
                        static void value(
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids)
                        {
                            info.limits->lock(lm_read_only);
                            data.u->lock(lm_read_and_write);
                            data.v->lock(lm_read_and_write);
                            data.h->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_only);
                            solids.solid_flags->lock(lm_read_only);
                            solids.solid_old_flags->lock(lm_read_and_write);

                            data.f_0->lock(lm_write_only);
                            data.f_1->lock(lm_write_only);
                            data.f_2->lock(lm_write_only);
                            data.f_3->lock(lm_write_only);
                            data.f_4->lock(lm_write_only);
                            data.f_5->lock(lm_write_only);
                            data.f_6->lock(lm_write_only);
                            data.f_7->lock(lm_write_only);
                            data.f_8->lock(lm_write_only);
                            data.f_eq_0->lock(lm_write_only);
                            data.f_eq_1->lock(lm_write_only);
                            data.f_eq_2->lock(lm_write_only);
                            data.f_eq_3->lock(lm_write_only);
                            data.f_eq_4->lock(lm_write_only);
                            data.f_eq_5->lock(lm_write_only);
                            data.f_eq_6->lock(lm_write_only);
                            data.f_eq_7->lock(lm_write_only);
                            data.f_eq_8->lock(lm_write_only);
                            data.f_temp_0->lock(lm_write_only);
                            data.f_temp_1->lock(lm_write_only);
                            data.f_temp_2->lock(lm_write_only);
                            data.f_temp_3->lock(lm_write_only);
                            data.f_temp_4->lock(lm_write_only);
                            data.f_temp_5->lock(lm_write_only);
                            data.f_temp_6->lock(lm_write_only);
                            data.f_temp_7->lock(lm_write_only);
                            data.f_temp_8->lock(lm_write_only);

                            info.cuda_dir_4->lock(lm_read_only);

                            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                            {
                                bool stf((*solids.solid_old_flags)[i] & !(*solids.solid_flags)[i]);
                                bool fts(!(*solids.solid_old_flags)[i] & (*solids.solid_flags)[i]);
                                (*solids.solid_old_flags)[i] = (*solids.solid_flags)[i];

                                (*data.h)[i] = stf ?
                                                    _extrapolation(info, data, solids, i) :(fts ? DT_(0) : (*data.h)[i]);

                                (*data.u)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_u  : (stf ?
                                                            _interpolation_u(info, data, solids, i) :(fts ? DT_(0) : (*data.u)[i]));

                                (*data.v)[i] = ((*solids.boundary_flags)[i] & !stf) ?
                                                    solids.current_v  : (stf ?
                                                            _interpolation_v(info, data, solids, i) :(fts ? DT_(0) : (*data.v)[i]));


                                (*data.f_0)[i] = fts ? DT_(0) : (*data.f_0)[i];
                                (*data.f_1)[i] = fts ? DT_(0) : (*data.f_1)[i];
                                (*data.f_2)[i] = fts ? DT_(0) : (*data.f_2)[i];
                                (*data.f_3)[i] = fts ? DT_(0) : (*data.f_3)[i];
                                (*data.f_4)[i] = fts ? DT_(0) : (*data.f_4)[i];
                                (*data.f_5)[i] = fts ? DT_(0) : (*data.f_5)[i];
                                (*data.f_6)[i] = fts ? DT_(0) : (*data.f_6)[i];
                                (*data.f_7)[i] = fts ? DT_(0) : (*data.f_7)[i];
                                (*data.f_8)[i] = fts ? DT_(0) : (*data.f_8)[i];

                                (*data.f_eq_0)[i] = fts ? DT_(0) : (*data.f_eq_0)[i];
                                (*data.f_eq_1)[i] = fts ? DT_(0) : (*data.f_eq_1)[i];
                                (*data.f_eq_2)[i] = fts ? DT_(0) : (*data.f_eq_2)[i];
                                (*data.f_eq_3)[i] = fts ? DT_(0) : (*data.f_eq_3)[i];
                                (*data.f_eq_4)[i] = fts ? DT_(0) : (*data.f_eq_4)[i];
                                (*data.f_eq_5)[i] = fts ? DT_(0) : (*data.f_eq_5)[i];
                                (*data.f_eq_6)[i] = fts ? DT_(0) : (*data.f_eq_6)[i];
                                (*data.f_eq_7)[i] = fts ? DT_(0) : (*data.f_eq_7)[i];
                                (*data.f_eq_8)[i] = fts ? DT_(0) : (*data.f_eq_8)[i];

                                (*data.f_temp_0)[i] = fts ? DT_(0) : (*data.f_temp_0)[i];
                                (*data.f_temp_1)[i] = fts ? DT_(0) : (*data.f_temp_1)[i];
                                (*data.f_temp_2)[i] = fts ? DT_(0) : (*data.f_temp_2)[i];
                                (*data.f_temp_3)[i] = fts ? DT_(0) : (*data.f_temp_3)[i];
                                (*data.f_temp_4)[i] = fts ? DT_(0) : (*data.f_temp_4)[i];
                                (*data.f_temp_5)[i] = fts ? DT_(0) : (*data.f_temp_5)[i];
                                (*data.f_temp_6)[i] = fts ? DT_(0) : (*data.f_temp_6)[i];
                                (*data.f_temp_7)[i] = fts ? DT_(0) : (*data.f_temp_7)[i];
                                (*data.f_temp_8)[i] = fts ? DT_(0) : (*data.f_temp_8)[i];
                            }

                            info.cuda_dir_4->unlock(lm_read_only);

                            solids.boundary_flags->unlock(lm_read_only);
                            solids.solid_flags->unlock(lm_read_only);
                            solids.solid_old_flags->unlock(lm_read_and_write);
                            data.h->unlock(lm_read_and_write);
                            data.v->unlock(lm_read_and_write);
                            data.u->unlock(lm_read_and_write);
                            info.limits->unlock(lm_read_only);

                            data.f_0->unlock(lm_write_only);
                            data.f_1->unlock(lm_write_only);
                            data.f_2->unlock(lm_write_only);
                            data.f_3->unlock(lm_write_only);
                            data.f_4->unlock(lm_write_only);
                            data.f_5->unlock(lm_write_only);
                            data.f_6->unlock(lm_write_only);
                            data.f_7->unlock(lm_write_only);
                            data.f_8->unlock(lm_write_only);
                            data.f_eq_0->unlock(lm_write_only);
                            data.f_eq_1->unlock(lm_write_only);
                            data.f_eq_2->unlock(lm_write_only);
                            data.f_eq_3->unlock(lm_write_only);
                            data.f_eq_4->unlock(lm_write_only);
                            data.f_eq_5->unlock(lm_write_only);
                            data.f_eq_6->unlock(lm_write_only);
                            data.f_eq_7->unlock(lm_write_only);
                            data.f_eq_8->unlock(lm_write_only);
                            data.f_temp_0->unlock(lm_write_only);
                            data.f_temp_1->unlock(lm_write_only);
                            data.f_temp_2->unlock(lm_write_only);
                            data.f_temp_3->unlock(lm_write_only);
                            data.f_temp_4->unlock(lm_write_only);
                            data.f_temp_5->unlock(lm_write_only);
                            data.f_temp_6->unlock(lm_write_only);
                            data.f_temp_7->unlock(lm_write_only);
                            data.f_temp_8->unlock(lm_write_only);
                        }
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_2>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_3>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_4>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_5>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_6>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_7>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
        template <>
            struct BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_8>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids);
            };
    }
}
#endif
