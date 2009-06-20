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

#ifndef LBM_GUARD_COLLIDE_STREAM_FSI_HH
#define LBM_GUARD_COLLIDE_STREAM_FSI_HH 1


#include <honei/lbm/tags.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <honei/util/benchmark_info.hh>
#include <cmath>

using namespace honei::lbm;

namespace honei
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
                          Grid<D2Q9, DT1_> & grid,
                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                          PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                          PackedSolidData<lbm_lattice_types::D2Q9, DT1_> & solids)
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

            //Determine cells with possibly needed backward-streaming
            std::vector<unsigned long> to_stream;
            for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
            {
                if((*solids.line_flags)[i] == true)
                {
                    to_stream.push_back(i);
                }
            }

            //Perform backward-streaming in all directions:
            for (unsigned long i(0) ; i < to_stream.size() ; ++i)
            {
                if((*data.f_temp_1)[to_stream[i]] != 0)
                {
                    (*data.f_temp_5)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i], solids.lines_inverse_j[i] - 1)] = (*data.f_temp_1)[to_stream[i]];
                    (*data.f_temp_1)[to_stream[i]] = DT1_(0);
                }

                if((*data.f_temp_2)[to_stream[i]] != 0)
                {
                    (*data.f_temp_6)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] - 1, solids.lines_inverse_j[i] - 1)] = (*data.f_temp_2)[to_stream[i]];
                    (*data.f_temp_2)[to_stream[i]] = DT1_(0);
                }

                if((*data.f_temp_3)[to_stream[i]] != 0)
                {
                    (*data.f_temp_7)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] - 1, solids.lines_inverse_j[i])] = (*data.f_temp_3)[to_stream[i]];
                    (*data.f_temp_3)[to_stream[i]] = DT1_(0);
                }


                if((*data.f_temp_4)[to_stream[i]] != 0)
                {
                    (*data.f_temp_8)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] - 1, solids.lines_inverse_j[i] + 1)] = (*data.f_temp_4)[to_stream[i]];
                    (*data.f_temp_4)[to_stream[i]] = DT1_(0);
                }

                if((*data.f_temp_5)[to_stream[i]] != 0)
                {
                    (*data.f_temp_1)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i], solids.lines_inverse_j[i] + 1)] = (*data.f_temp_5)[to_stream[i]];
                    (*data.f_temp_5)[to_stream[i]] = DT1_(0);
                }


                if((*data.f_temp_6)[to_stream[i]] != 0)
                {
                    (*data.f_temp_2)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] + 1, solids.lines_inverse_j[i] + 1)] = (*data.f_temp_6)[to_stream[i]];
                    (*data.f_temp_6)[to_stream[i]] = DT1_(0);
                }

                if((*data.f_temp_7)[to_stream[i]] != 0)
                {
                    (*data.f_temp_3)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] + 1, solids.lines_inverse_j[i])] = (*data.f_temp_7)[to_stream[i]];
                    (*data.f_temp_7)[to_stream[i]] = DT1_(0);
                }

                if((*data.f_temp_8)[to_stream[i]] != 0)
                {
                    (*data.f_temp_4)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT1_>::
                        h_index(grid, solids.lines_inverse_i[i] + 1, solids.lines_inverse_j[i] - 1)] = (*data.f_temp_8)[to_stream[i]];
                    (*data.f_temp_8)[to_stream[i]] = DT1_(0);
                }
            }

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
        }

    };
}
#endif
