/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008, 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LBM_GUARD_COLLIDE_STREAM_GRID_HH
#define LBM_GUARD_COLLIDE_STREAM_GRID_HH 1


/**
 * \file
 * Implementation of collision and streaming modules used by  LBM - (SWE, NavSto) solvers (using PackedGrid).
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <honei/util/benchmark_info.hh>
#include <cmath>

using namespace honei::lbm;

namespace honei
{
    template <typename Tag_, typename BoundaryType_, typename LatticeType_>
    struct CollideStreamGrid
    {
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <>
    struct CollideStreamGrid<tags::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        /**
         * \name Collision and Streaming for direction 1..
         *
         * \brief Solves the LB equation.
         *
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(
                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                          PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming:");

            info.limits->lock(lm_read_only);
            info.dir_1->lock(lm_read_only);
            info.dir_2->lock(lm_read_only);
            info.dir_3->lock(lm_read_only);
            info.dir_4->lock(lm_read_only);
            info.dir_5->lock(lm_read_only);
            info.dir_6->lock(lm_read_only);
            info.dir_7->lock(lm_read_only);
            info.dir_8->lock(lm_read_only);
            info.dir_index_1->lock(lm_read_only);
            info.dir_index_2->lock(lm_read_only);
            info.dir_index_3->lock(lm_read_only);
            info.dir_index_4->lock(lm_read_only);
            info.dir_index_5->lock(lm_read_only);
            info.dir_index_6->lock(lm_read_only);
            info.dir_index_7->lock(lm_read_only);
            info.dir_index_8->lock(lm_read_only);

            data.f_eq_0->lock(lm_read_only);
            data.f_eq_1->lock(lm_read_only);
            data.f_eq_2->lock(lm_read_only);
            data.f_eq_3->lock(lm_read_only);
            data.f_eq_4->lock(lm_read_only);
            data.f_eq_5->lock(lm_read_only);
            data.f_eq_6->lock(lm_read_only);
            data.f_eq_7->lock(lm_read_only);
            data.f_eq_8->lock(lm_read_only);
            data.f_0->lock(lm_read_only);
            data.f_1->lock(lm_read_only);
            data.f_2->lock(lm_read_only);
            data.f_3->lock(lm_read_only);
            data.f_4->lock(lm_read_only);
            data.f_5->lock(lm_read_only);
            data.f_6->lock(lm_read_only);
            data.f_7->lock(lm_read_only);
            data.f_8->lock(lm_read_only);

            data.f_temp_0->lock(lm_write_only);
            data.f_temp_1->lock(lm_write_only);
            data.f_temp_2->lock(lm_write_only);
            data.f_temp_3->lock(lm_write_only);
            data.f_temp_4->lock(lm_write_only);
            data.f_temp_5->lock(lm_write_only);
            data.f_temp_6->lock(lm_write_only);
            data.f_temp_7->lock(lm_write_only);
            data.f_temp_8->lock(lm_write_only);

            for (unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() - 1] ; ++i)
            {
                (*data.f_temp_0)[i] = (*data.f_0)[i] - ((*data.f_0)[i] - (*data.f_eq_0)[i])/tau;
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_1)[(*info.dir_1)[half] + offset] = (*data.f_1)[i] - ((*data.f_1)[i] - (*data.f_eq_1)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_2)[(*info.dir_2)[half] + offset] = (*data.f_2)[i] - ((*data.f_2)[i] - (*data.f_eq_2)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_3)[(*info.dir_3)[half] + offset] = (*data.f_3)[i] - ((*data.f_3)[i] - (*data.f_eq_3)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_4)[(*info.dir_4)[half] + offset] = (*data.f_4)[i] - ((*data.f_4)[i] - (*data.f_eq_4)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_5)[(*info.dir_5)[half] + offset] = (*data.f_5)[i] - ((*data.f_5)[i] - (*data.f_eq_5)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_6)[(*info.dir_6)[half] + offset] = (*data.f_6)[i] - ((*data.f_6)[i] - (*data.f_eq_6)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_7)[(*info.dir_7)[half] + offset] = (*data.f_7)[i] - ((*data.f_7)[i] - (*data.f_eq_7)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_8)[(*info.dir_8)[half] + offset] = (*data.f_8)[i] - ((*data.f_8)[i] - (*data.f_eq_8)[i])/tau;
                }
            }

            info.limits->unlock(lm_read_only);
            info.dir_1->unlock(lm_read_only);
            info.dir_2->unlock(lm_read_only);
            info.dir_3->unlock(lm_read_only);
            info.dir_4->unlock(lm_read_only);
            info.dir_5->unlock(lm_read_only);
            info.dir_6->unlock(lm_read_only);
            info.dir_7->unlock(lm_read_only);
            info.dir_8->unlock(lm_read_only);
            info.dir_index_1->unlock(lm_read_only);
            info.dir_index_2->unlock(lm_read_only);
            info.dir_index_3->unlock(lm_read_only);
            info.dir_index_4->unlock(lm_read_only);
            info.dir_index_5->unlock(lm_read_only);
            info.dir_index_6->unlock(lm_read_only);
            info.dir_index_7->unlock(lm_read_only);
            info.dir_index_8->unlock(lm_read_only);

            data.f_eq_0->unlock(lm_read_only);
            data.f_eq_1->unlock(lm_read_only);
            data.f_eq_2->unlock(lm_read_only);
            data.f_eq_3->unlock(lm_read_only);
            data.f_eq_4->unlock(lm_read_only);
            data.f_eq_5->unlock(lm_read_only);
            data.f_eq_6->unlock(lm_read_only);
            data.f_eq_7->unlock(lm_read_only);
            data.f_eq_8->unlock(lm_read_only);
            data.f_0->unlock(lm_read_only);
            data.f_1->unlock(lm_read_only);
            data.f_2->unlock(lm_read_only);
            data.f_3->unlock(lm_read_only);
            data.f_4->unlock(lm_read_only);
            data.f_5->unlock(lm_read_only);
            data.f_6->unlock(lm_read_only);
            data.f_7->unlock(lm_read_only);
            data.f_8->unlock(lm_read_only);

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

        template<typename DT1_>
            static inline BenchmarkInfo get_benchmark_info(HONEI_UNUSED PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
            {
                BenchmarkInfo result;
                result.flops = data->h->size() * 9 * 3;
                result.load = data->h->size() * 9 * 4 * sizeof(DT1_);
                result.store = data->h->size() * 9 * sizeof(DT1_);
                result.size.push_back(data->h->size());
                return result;
            }
    };

    template <>
    struct CollideStreamGrid<tags::CPU::Generic, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        template <typename DT1_, typename DT2_>
        static void value(
                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                          PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming:");

            info.limits->lock(lm_read_only);
            info.dir_1->lock(lm_read_only);
            info.dir_2->lock(lm_read_only);
            info.dir_3->lock(lm_read_only);
            info.dir_4->lock(lm_read_only);
            info.dir_5->lock(lm_read_only);
            info.dir_6->lock(lm_read_only);
            info.dir_7->lock(lm_read_only);
            info.dir_8->lock(lm_read_only);
            info.dir_index_1->lock(lm_read_only);
            info.dir_index_2->lock(lm_read_only);
            info.dir_index_3->lock(lm_read_only);
            info.dir_index_4->lock(lm_read_only);
            info.dir_index_5->lock(lm_read_only);
            info.dir_index_6->lock(lm_read_only);
            info.dir_index_7->lock(lm_read_only);
            info.dir_index_8->lock(lm_read_only);

            data.f_eq_0->lock(lm_read_only);
            data.f_eq_1->lock(lm_read_only);
            data.f_eq_2->lock(lm_read_only);
            data.f_eq_3->lock(lm_read_only);
            data.f_eq_4->lock(lm_read_only);
            data.f_eq_5->lock(lm_read_only);
            data.f_eq_6->lock(lm_read_only);
            data.f_eq_7->lock(lm_read_only);
            data.f_eq_8->lock(lm_read_only);
            data.f_0->lock(lm_read_only);
            data.f_1->lock(lm_read_only);
            data.f_2->lock(lm_read_only);
            data.f_3->lock(lm_read_only);
            data.f_4->lock(lm_read_only);
            data.f_5->lock(lm_read_only);
            data.f_6->lock(lm_read_only);
            data.f_7->lock(lm_read_only);
            data.f_8->lock(lm_read_only);

            data.f_temp_0->lock(lm_write_only);
            data.f_temp_1->lock(lm_write_only);
            data.f_temp_2->lock(lm_write_only);
            data.f_temp_3->lock(lm_write_only);
            data.f_temp_4->lock(lm_write_only);
            data.f_temp_5->lock(lm_write_only);
            data.f_temp_6->lock(lm_write_only);
            data.f_temp_7->lock(lm_write_only);
            data.f_temp_8->lock(lm_write_only);

            const unsigned long * const limits(info.limits->elements());
            const unsigned long * const dir_1(info.dir_1->elements());
            const unsigned long * const dir_2(info.dir_2->elements());
            const unsigned long * const dir_3(info.dir_3->elements());
            const unsigned long * const dir_4(info.dir_4->elements());
            const unsigned long * const dir_5(info.dir_5->elements());
            const unsigned long * const dir_6(info.dir_6->elements());
            const unsigned long * const dir_7(info.dir_7->elements());
            const unsigned long * const dir_8(info.dir_8->elements());
            const unsigned long * const dir_index_1(info.dir_index_1->elements());
            const unsigned long * const dir_index_2(info.dir_index_2->elements());
            const unsigned long * const dir_index_3(info.dir_index_3->elements());
            const unsigned long * const dir_index_4(info.dir_index_4->elements());
            const unsigned long * const dir_index_5(info.dir_index_5->elements());
            const unsigned long * const dir_index_6(info.dir_index_6->elements());
            const unsigned long * const dir_index_7(info.dir_index_7->elements());
            const unsigned long * const dir_index_8(info.dir_index_8->elements());

            const DT1_ * const f_eq_0(data.f_eq_0->elements());
            const DT1_ * const f_eq_1(data.f_eq_1->elements());
            const DT1_ * const f_eq_2(data.f_eq_2->elements());
            const DT1_ * const f_eq_3(data.f_eq_3->elements());
            const DT1_ * const f_eq_4(data.f_eq_4->elements());
            const DT1_ * const f_eq_5(data.f_eq_5->elements());
            const DT1_ * const f_eq_6(data.f_eq_6->elements());
            const DT1_ * const f_eq_7(data.f_eq_7->elements());
            const DT1_ * const f_eq_8(data.f_eq_8->elements());
            const DT1_ * const f_0(data.f_0->elements());
            const DT1_ * const f_1(data.f_1->elements());
            const DT1_ * const f_2(data.f_2->elements());
            const DT1_ * const f_3(data.f_3->elements());
            const DT1_ * const f_4(data.f_4->elements());
            const DT1_ * const f_5(data.f_5->elements());
            const DT1_ * const f_6(data.f_6->elements());
            const DT1_ * const f_7(data.f_7->elements());
            const DT1_ * const f_8(data.f_8->elements());

            DT1_ * f_temp_0(data.f_temp_0->elements());
            DT1_ * f_temp_1(data.f_temp_1->elements());
            DT1_ * f_temp_2(data.f_temp_2->elements());
            DT1_ * f_temp_3(data.f_temp_3->elements());
            DT1_ * f_temp_4(data.f_temp_4->elements());
            DT1_ * f_temp_5(data.f_temp_5->elements());
            DT1_ * f_temp_6(data.f_temp_6->elements());
            DT1_ * f_temp_7(data.f_temp_7->elements());
            DT1_ * f_temp_8(data.f_temp_8->elements());

            for (unsigned long i((limits)[0]) ; i < (limits)[info.limits->size() - 1] ; ++i)
            {
                (f_temp_0)[i] = (f_0)[i] - ((f_0)[i] - (f_eq_0)[i])/tau;
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_1[begin + 1]);
                for (unsigned long i((dir_index_1)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_1)[(dir_1)[half] + offset] = (f_1)[i] - ((f_1)[i] - (f_eq_1)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_2[begin + 1]);
                for (unsigned long i((dir_index_2)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_2)[(dir_2)[half] + offset] = (f_2)[i] - ((f_2)[i] - (f_eq_2)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_3[begin + 1]);
                for (unsigned long i((dir_index_3)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_3)[(dir_3)[half] + offset] = (f_3)[i] - ((f_3)[i] - (f_eq_3)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_4[begin + 1]);
                for (unsigned long i((dir_index_4)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_4)[(dir_4)[half] + offset] = (f_4)[i] - ((f_4)[i] - (f_eq_4)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_5[begin + 1]);
                for (unsigned long i((dir_index_5)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_5)[(dir_5)[half] + offset] = (f_5)[i] - ((f_5)[i] - (f_eq_5)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_6[begin + 1]);
                for (unsigned long i((dir_index_6)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_6)[(dir_6)[half] + offset] = (f_6)[i] - ((f_6)[i] - (f_eq_6)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_7[begin + 1]);
                for (unsigned long i((dir_index_7)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_7)[(dir_7)[half] + offset] = (f_7)[i] - ((f_7)[i] - (f_eq_7)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
            {
                const unsigned long end(dir_index_8[begin + 1]);
                for (unsigned long i((dir_index_8)[begin]), offset(0) ; i < end ; ++i, ++offset)
                {
                    (f_temp_8)[(dir_8)[half] + offset] = (f_8)[i] - ((f_8)[i] - (f_eq_8)[i])/tau;
                }
            }

            info.limits->unlock(lm_read_only);
            info.dir_1->unlock(lm_read_only);
            info.dir_2->unlock(lm_read_only);
            info.dir_3->unlock(lm_read_only);
            info.dir_4->unlock(lm_read_only);
            info.dir_5->unlock(lm_read_only);
            info.dir_6->unlock(lm_read_only);
            info.dir_7->unlock(lm_read_only);
            info.dir_8->unlock(lm_read_only);
            info.dir_index_1->unlock(lm_read_only);
            info.dir_index_2->unlock(lm_read_only);
            info.dir_index_3->unlock(lm_read_only);
            info.dir_index_4->unlock(lm_read_only);
            info.dir_index_5->unlock(lm_read_only);
            info.dir_index_6->unlock(lm_read_only);
            info.dir_index_7->unlock(lm_read_only);
            info.dir_index_8->unlock(lm_read_only);

            data.f_eq_0->unlock(lm_read_only);
            data.f_eq_1->unlock(lm_read_only);
            data.f_eq_2->unlock(lm_read_only);
            data.f_eq_3->unlock(lm_read_only);
            data.f_eq_4->unlock(lm_read_only);
            data.f_eq_5->unlock(lm_read_only);
            data.f_eq_6->unlock(lm_read_only);
            data.f_eq_7->unlock(lm_read_only);
            data.f_eq_8->unlock(lm_read_only);
            data.f_0->unlock(lm_read_only);
            data.f_1->unlock(lm_read_only);
            data.f_2->unlock(lm_read_only);
            data.f_3->unlock(lm_read_only);
            data.f_4->unlock(lm_read_only);
            data.f_5->unlock(lm_read_only);
            data.f_6->unlock(lm_read_only);
            data.f_7->unlock(lm_read_only);
            data.f_8->unlock(lm_read_only);

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
    struct CollideStreamGrid<tags::GPU::CUDA, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                float tau);
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, double> & data,
                double tau);
    };

    template <>
    struct CollideStreamGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        template <typename DT_>
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                DT_ tau);
    };

    template <>
    struct CollideStreamGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        template <typename DT_>
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                DT_ tau);
    };

    template <>
    struct CollideStreamGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        template <typename DT1_>
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                DT1_ tau);
    };

    template <>
    struct CollideStreamGrid<tags::Cell, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                float tau);
    };

    template <>
    struct CollideStreamGrid<tags::CPU::Itanium, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        template <typename DT1_>
        static void value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                DT1_ tau);
    };
}
#endif
