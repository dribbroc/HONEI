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
#ifndef LBM_GUARD_COLLIDE_STREAM_FSI_HH
#define LBM_GUARD_COLLIDE_STREAM_FSI_HH 1


#include <honei/lbm/tags.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>

using namespace honei::lbm;

namespace honei
{
    namespace lbm
    {
        template <typename Tag_, typename BoundaryType_, typename LatticeType_>
            struct CollideStreamFSI
            {
            };

        /**
         * \brief Collision and streaming module for LABSWE.
         *
         * \ingroup grplbmoperations
         */
        template <>
            struct CollideStreamFSI<tags::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
            {
                template <typename DT1_>
                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, DT1_> & solids,
                            DT1_ d_x,
                            DT1_ d_y)
                    {
                        CONTEXT("When performing collision and streaming FSI correction:");

                        info.limits->lock(lm_read_only);
                        data.f_temp_0->lock(lm_read_and_write);
                        data.f_temp_1->lock(lm_read_and_write);
                        data.f_temp_2->lock(lm_read_and_write);
                        data.f_temp_3->lock(lm_read_and_write);
                        data.f_temp_4->lock(lm_read_and_write);
                        data.f_temp_5->lock(lm_read_and_write);
                        data.f_temp_6->lock(lm_read_and_write);
                        data.f_temp_7->lock(lm_read_and_write);
                        data.f_temp_8->lock(lm_read_and_write);

                        info.cuda_dir_1->lock(lm_read_only);
                        info.cuda_dir_2->lock(lm_read_only);
                        info.cuda_dir_3->lock(lm_read_only);
                        info.cuda_dir_4->lock(lm_read_only);
                        info.cuda_dir_5->lock(lm_read_only);
                        info.cuda_dir_6->lock(lm_read_only);
                        info.cuda_dir_7->lock(lm_read_only);
                        info.cuda_dir_8->lock(lm_read_only);

                        solids.solid_flags->lock(lm_read_only);

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        solids.f_mea_1->lock(lm_read_and_write);
                        solids.f_mea_2->lock(lm_read_and_write);
                        solids.f_mea_3->lock(lm_read_and_write);
                        solids.f_mea_4->lock(lm_read_and_write);
                        solids.f_mea_5->lock(lm_read_and_write);
                        solids.f_mea_6->lock(lm_read_and_write);
                        solids.f_mea_7->lock(lm_read_and_write);
                        solids.f_mea_8->lock(lm_read_and_write);

                        //Perform backward-streaming in all directions:
                        for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
                        {
                            bool valid(((*data.f_temp_1)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            unsigned long prev_index(valid ? (*info.cuda_dir_5)[i] : i);
                            DT1_ mb_term(DT1_(6. * 1./9.) * (d_x * solids.current_u * (*data.distribution_x)[5] +
                                                             d_y * solids.current_v * (*data.distribution_y)[5]));
                            (*data.f_temp_5)[prev_index] = valid ? (*data.f_temp_1)[prev_index] + mb_term : (*data.f_temp_5)[prev_index];
                            (*solids.f_mea_5)[i] = valid ? (((*data.distribution_x)[5] + (*data.distribution_y)[5]) *
                                                                ((*data.f_temp_5)[prev_index] + (*data.f_temp_1)[i])) : (*solids.f_mea_5)[prev_index];
                            (*data.f_temp_1)[i] = valid ? DT1_(0) : (*data.f_temp_1)[i];

                            valid = (((*data.f_temp_2)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_6)[i] : i);
                            mb_term = (DT1_(6. * 1./36.) * (d_x * solids.current_u * (*data.distribution_x)[6] +
                                                             d_y * solids.current_v * (*data.distribution_y)[6]));
                            (*data.f_temp_6)[prev_index] = valid ? (*data.f_temp_2)[prev_index] + mb_term : (*data.f_temp_6)[prev_index];
                            (*solids.f_mea_6)[i] = valid ? (((*data.distribution_x)[6] + (*data.distribution_y)[6]) *
                                                                ((*data.f_temp_6)[prev_index] + (*data.f_temp_2)[i])) : (*solids.f_mea_6)[prev_index];
                            (*data.f_temp_2)[i] = valid ? DT1_(0) : (*data.f_temp_2)[i];

                            valid = (((*data.f_temp_3)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_7)[i] : i);
                            mb_term = (DT1_(6. * 1./9.) * (d_x * solids.current_u * (*data.distribution_x)[7] +
                                                             d_y * solids.current_v * (*data.distribution_y)[7]));
                            (*data.f_temp_7)[prev_index] = valid ? (*data.f_temp_3)[prev_index] + mb_term : (*data.f_temp_7)[prev_index];
                            (*solids.f_mea_7)[i] = valid ? (((*data.distribution_x)[7] + (*data.distribution_y)[7]) *
                                                                ((*data.f_temp_7)[prev_index] + (*data.f_temp_3)[i])) : (*solids.f_mea_7)[prev_index];
                            (*data.f_temp_3)[i] = valid ? DT1_(0) : (*data.f_temp_3)[i];

                            valid = (((*data.f_temp_4)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_8)[i] : i);
                            mb_term = (DT1_(6. * 1./36.) * (d_x * solids.current_u * (*data.distribution_x)[8] +
                                                             d_y * solids.current_v * (*data.distribution_y)[8]));
                            (*data.f_temp_8)[prev_index] = valid ? (*data.f_temp_4)[prev_index] + mb_term : (*data.f_temp_8)[prev_index];
                            (*solids.f_mea_8)[i] = valid ? (((*data.distribution_x)[8] + (*data.distribution_y)[8]) *
                                                                ((*data.f_temp_8)[prev_index] + (*data.f_temp_4)[i])) : (*solids.f_mea_8)[prev_index];
                            (*data.f_temp_4)[i] = valid ? DT1_(0) : (*data.f_temp_4)[i];

                            valid = (((*data.f_temp_5)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_1)[i] : i);
                            mb_term = (DT1_(6. * 1./9.) * (d_x * solids.current_u * (*data.distribution_x)[1] +
                                                             d_y * solids.current_v * (*data.distribution_y)[1]));
                            (*data.f_temp_1)[prev_index] = valid ? (*data.f_temp_5)[prev_index] + mb_term : (*data.f_temp_1)[prev_index];
                            (*solids.f_mea_1)[i] = valid ? (((*data.distribution_x)[1] + (*data.distribution_y)[1]) *
                                                                ((*data.f_temp_1)[prev_index] + (*data.f_temp_5)[i])) : (*solids.f_mea_1)[prev_index];
                            (*data.f_temp_5)[i] = valid ? DT1_(0) : (*data.f_temp_5)[i];

                            valid = (((*data.f_temp_6)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_2)[i] : i);
                            mb_term = (DT1_(6. * 1./36.) * (d_x * solids.current_u * (*data.distribution_x)[2] +
                                                             d_y * solids.current_v * (*data.distribution_y)[2]));
                            (*data.f_temp_2)[prev_index] = valid ? (*data.f_temp_6)[prev_index] + mb_term : (*data.f_temp_2)[prev_index];
                            (*solids.f_mea_2)[i] = valid ? (((*data.distribution_x)[2] + (*data.distribution_y)[2]) *
                                                                ((*data.f_temp_2)[prev_index] + (*data.f_temp_6)[i])) : (*solids.f_mea_2)[prev_index];
                            (*data.f_temp_6)[i] = valid ? DT1_(0) : (*data.f_temp_6)[i];

                            valid = (((*data.f_temp_7)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_3)[i] : i);
                            mb_term = (DT1_(6. * 1./9.) * (d_x * solids.current_u * (*data.distribution_x)[3] +
                                                             d_y * solids.current_v * (*data.distribution_y)[3]));
                            (*data.f_temp_3)[prev_index] = valid ? (*data.f_temp_7)[prev_index] + mb_term : (*data.f_temp_3)[prev_index];
                            (*solids.f_mea_3)[i] = valid ? (((*data.distribution_x)[3] + (*data.distribution_y)[3]) *
                                                                ((*data.f_temp_3)[prev_index] + (*data.f_temp_7)[i])) : (*solids.f_mea_3)[prev_index];
                            (*data.f_temp_7)[i] = valid ? DT1_(0) : (*data.f_temp_7)[i];

                            valid = (((*data.f_temp_8)[i] != DT1_(0) && (*solids.solid_flags)[i]));
                            prev_index = (valid ? (*info.cuda_dir_4)[i] : i);
                            mb_term = (DT1_(6. * 1./36.) * (d_x * solids.current_u * (*data.distribution_x)[4] +
                                                             d_y * solids.current_v * (*data.distribution_y)[4]));
                            (*data.f_temp_4)[prev_index] = valid ? (*data.f_temp_8)[prev_index] + mb_term : (*data.f_temp_4)[prev_index];
                            (*solids.f_mea_4)[i] = valid ? (((*data.distribution_x)[4] + (*data.distribution_y)[4]) *
                                                                ((*data.f_temp_4)[prev_index] + (*data.f_temp_8)[i])) : (*solids.f_mea_4)[prev_index];
                            (*data.f_temp_8)[i] = valid ? DT1_(0) : (*data.f_temp_8)[i];
                        }

                        solids.solid_flags->unlock(lm_read_only);

                        info.cuda_dir_1->unlock(lm_read_only);
                        info.cuda_dir_2->unlock(lm_read_only);
                        info.cuda_dir_3->unlock(lm_read_only);
                        info.cuda_dir_4->unlock(lm_read_only);
                        info.cuda_dir_5->unlock(lm_read_only);
                        info.cuda_dir_6->unlock(lm_read_only);
                        info.cuda_dir_7->unlock(lm_read_only);
                        info.cuda_dir_8->unlock(lm_read_only);

                        info.limits->unlock(lm_read_only);
                        data.f_temp_0->unlock(lm_read_and_write);
                        data.f_temp_1->unlock(lm_read_and_write);
                        data.f_temp_2->unlock(lm_read_and_write);
                        data.f_temp_3->unlock(lm_read_and_write);
                        data.f_temp_4->unlock(lm_read_and_write);
                        data.f_temp_5->unlock(lm_read_and_write);
                        data.f_temp_6->unlock(lm_read_and_write);
                        data.f_temp_7->unlock(lm_read_and_write);
                        data.f_temp_8->unlock(lm_read_and_write);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        solids.f_mea_1->unlock(lm_read_and_write);
                        solids.f_mea_2->unlock(lm_read_and_write);
                        solids.f_mea_3->unlock(lm_read_and_write);
                        solids.f_mea_4->unlock(lm_read_and_write);
                        solids.f_mea_5->unlock(lm_read_and_write);
                        solids.f_mea_6->unlock(lm_read_and_write);
                        solids.f_mea_7->unlock(lm_read_and_write);
                        solids.f_mea_8->unlock(lm_read_and_write);
                    }

            };

        template <>
            struct CollideStreamFSI<tags::GPU::CUDA, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
            {

                    static void value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids,
                            float d_x,
                            float d_y);
            };


    }
}
#endif
