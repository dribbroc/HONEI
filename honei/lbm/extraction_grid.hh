/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
#ifndef LBM_GUARD_EXTRACTION_HH
#define LBM_GUARD_EXTRACTION_HH 1


/**
 * \file
 * Implementations of quantity extraction modules used by LBM - (SWE) solvers using a PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/lbm_limiter.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/attributes.hh>

using namespace lbm;

namespace honei
{
    template<typename Tag_, typename LbmMode_>
        struct ExtractionGrid
        {
        };

    template<>
        struct ExtractionGrid<tags::CPU, lbm_modes::DRY>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ epsilon)
                    {
                        CONTEXT("When extracting h, u and v:");

                        //set f to t_temp
                        DenseVector<DT_> * swap;
                        swap = data.f_0;
                        data.f_0 = data.f_temp_0;
                        data.f_temp_0 = swap;
                        swap = data.f_1;
                        data.f_1 = data.f_temp_1;
                        data.f_temp_1 = swap;
                        swap = data.f_2;
                        data.f_2 = data.f_temp_2;
                        data.f_temp_2 = swap;
                        swap = data.f_3;
                        data.f_3 = data.f_temp_3;
                        data.f_temp_3 = swap;
                        swap = data.f_4;
                        data.f_4 = data.f_temp_4;
                        data.f_temp_4 = swap;
                        swap = data.f_5;
                        data.f_5 = data.f_temp_5;
                        data.f_temp_5 = swap;
                        swap = data.f_6;
                        data.f_6 = data.f_temp_6;
                        data.f_temp_6 = swap;
                        swap = data.f_7;
                        data.f_7 = data.f_temp_7;
                        data.f_temp_7 = swap;
                        swap = data.f_8;
                        data.f_8 = data.f_temp_8;
                        data.f_temp_8 = swap;

                        info.limits->lock(lm_read_only);

                        data.f_0->lock(lm_read_and_write);
                        data.f_1->lock(lm_read_and_write);
                        data.f_2->lock(lm_read_and_write);
                        data.f_3->lock(lm_read_and_write);
                        data.f_4->lock(lm_read_and_write);
                        data.f_5->lock(lm_read_and_write);
                        data.f_6->lock(lm_read_and_write);
                        data.f_7->lock(lm_read_and_write);
                        data.f_8->lock(lm_read_and_write);

                        data.h->lock(lm_write_only);

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        data.u->lock(lm_write_only);
                        data.v->lock(lm_write_only);

                        DT_ lax_upper(epsilon /*10e-5 std::numeric_limits<DT_>::epsilon() * 10e7*/);
                        DT_ lax_lower(-lax_upper);

                        for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                        {

                            //accumulate
                            (*data.h)[i] = (*data.f_0)[i] +
                                (*data.f_1)[i] +
                                (*data.f_2)[i] +
                                (*data.f_3)[i] +
                                (*data.f_4)[i] +
                                (*data.f_5)[i] +
                                (*data.f_6)[i] +
                                (*data.f_7)[i] +
                                (*data.f_8)[i];

                            if((*data.h)[i] < lax_lower || (*data.h)[i] > lax_upper)
                            {
                                (*data.u)[i] = ((*data.distribution_x)[0] * (*data.f_0)[i] +
                                        (*data.distribution_x)[1] * (*data.f_1)[i] +
                                        (*data.distribution_x)[2] * (*data.f_2)[i] +
                                        (*data.distribution_x)[3] * (*data.f_3)[i] +
                                        (*data.distribution_x)[4] * (*data.f_4)[i] +
                                        (*data.distribution_x)[5] * (*data.f_5)[i] +
                                        (*data.distribution_x)[6] * (*data.f_6)[i] +
                                        (*data.distribution_x)[7] * (*data.f_7)[i] +
                                        (*data.distribution_x)[8] * (*data.f_8)[i]) / (*data.h)[i];

                                (*data.v)[i] = ((*data.distribution_y)[0] * (*data.f_0)[i] +
                                        (*data.distribution_y)[1] * (*data.f_1)[i] +
                                        (*data.distribution_y)[2] * (*data.f_2)[i] +
                                        (*data.distribution_y)[3] * (*data.f_3)[i] +
                                        (*data.distribution_y)[4] * (*data.f_4)[i] +
                                        (*data.distribution_y)[5] * (*data.f_5)[i] +
                                        (*data.distribution_y)[6] * (*data.f_6)[i] +
                                        (*data.distribution_y)[7] * (*data.f_7)[i] +
                                        (*data.distribution_y)[8] * (*data.f_8)[i]) / (*data.h)[i];
                            }
                            else
                            {
                                //TODO: better heuristics for reset of h: if negative -> 0, if positive but too small -> epsilon
                                (*data.h)[i] = 0;

                                (*data.u)[i] = 0;
                                (*data.v)[i] = 0;
                            }
                            (*data.h)[i] = MinModLimiter<tags::CPU>::value((*data.h)[i]);


                        }

                        info.limits->unlock(lm_read_only);

                        data.f_0->unlock(lm_read_and_write);
                        data.f_1->unlock(lm_read_and_write);
                        data.f_2->unlock(lm_read_and_write);
                        data.f_3->unlock(lm_read_and_write);
                        data.f_4->unlock(lm_read_and_write);
                        data.f_5->unlock(lm_read_and_write);
                        data.f_6->unlock(lm_read_and_write);
                        data.f_7->unlock(lm_read_and_write);
                        data.f_8->unlock(lm_read_and_write);

                        data.h->unlock(lm_write_only);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        data.u->unlock(lm_write_only);
                        data.v->unlock(lm_write_only);
                    }

                template<typename DT1_>
                    static inline BenchmarkInfo get_benchmark_info(HONEI_UNUSED PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                    {
                        BenchmarkInfo result;
                        result.flops = data->h->size() * 44;
                        result.load = data->h->size() * (5 * 9 + 2) * sizeof(DT1_);
                        result.store = data->h->size() * 3 * sizeof(DT1_);
                        result.size.push_back(data->h->size());
                        return result;
                    }
        };

    template<>
        struct ExtractionGrid<tags::CPU::Generic, lbm_modes::DRY>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ epsilon)
                    {
                        CONTEXT("When extracting h, u and v:");

                        //set f to t_temp
                        DenseVector<DT_> * swap;
                        swap = data.f_0;
                        data.f_0 = data.f_temp_0;
                        data.f_temp_0 = swap;
                        swap = data.f_1;
                        data.f_1 = data.f_temp_1;
                        data.f_temp_1 = swap;
                        swap = data.f_2;
                        data.f_2 = data.f_temp_2;
                        data.f_temp_2 = swap;
                        swap = data.f_3;
                        data.f_3 = data.f_temp_3;
                        data.f_temp_3 = swap;
                        swap = data.f_4;
                        data.f_4 = data.f_temp_4;
                        data.f_temp_4 = swap;
                        swap = data.f_5;
                        data.f_5 = data.f_temp_5;
                        data.f_temp_5 = swap;
                        swap = data.f_6;
                        data.f_6 = data.f_temp_6;
                        data.f_temp_6 = swap;
                        swap = data.f_7;
                        data.f_7 = data.f_temp_7;
                        data.f_temp_7 = swap;
                        swap = data.f_8;
                        data.f_8 = data.f_temp_8;
                        data.f_temp_8 = swap;

                        info.limits->lock(lm_read_only);

                        data.f_0->lock(lm_read_and_write);
                        data.f_1->lock(lm_read_and_write);
                        data.f_2->lock(lm_read_and_write);
                        data.f_3->lock(lm_read_and_write);
                        data.f_4->lock(lm_read_and_write);
                        data.f_5->lock(lm_read_and_write);
                        data.f_6->lock(lm_read_and_write);
                        data.f_7->lock(lm_read_and_write);
                        data.f_8->lock(lm_read_and_write);

                        data.h->lock(lm_write_only);

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        data.u->lock(lm_write_only);
                        data.v->lock(lm_write_only);

                        const unsigned long * const limits(info.limits->elements());
                        DT_ * f_0(data.f_0->elements());
                        DT_ * f_1(data.f_1->elements());
                        DT_ * f_2(data.f_2->elements());
                        DT_ * f_3(data.f_3->elements());
                        DT_ * f_4(data.f_4->elements());
                        DT_ * f_5(data.f_5->elements());
                        DT_ * f_6(data.f_6->elements());
                        DT_ * f_7(data.f_7->elements());
                        DT_ * f_8(data.f_8->elements());

                        DT_ * h(data.h->elements());
                        DT_ * u(data.u->elements());
                        DT_ * v(data.v->elements());
                        const DT_ * const distribution_x(data.distribution_x->elements());
                        const DT_ * const distribution_y(data.distribution_y->elements());

                        DT_ lax_upper(epsilon /*10e-5 std::numeric_limits<DT_>::epsilon() * 10e7*/);
                        DT_ lax_lower(-lax_upper);

                        const unsigned long start(limits[0]);
                        const unsigned long end(limits[info.limits->size() - 1]);
                        for(unsigned long i(start); i < end; ++i)
                        {

                            //accumulate
                            (h)[i] = (f_0)[i] +
                                (f_1)[i] +
                                (f_2)[i] +
                                (f_3)[i] +
                                (f_4)[i] +
                                (f_5)[i] +
                                (f_6)[i] +
                                (f_7)[i] +
                                (f_8)[i];

                            if((h)[i] < lax_lower || (h)[i] > lax_upper)
                            {
                                (u)[i] = ((distribution_x)[0] * (f_0)[i] +
                                        (distribution_x)[1] * (f_1)[i] +
                                        (distribution_x)[2] * (f_2)[i] +
                                        (distribution_x)[3] * (f_3)[i] +
                                        (distribution_x)[4] * (f_4)[i] +
                                        (distribution_x)[5] * (f_5)[i] +
                                        (distribution_x)[6] * (f_6)[i] +
                                        (distribution_x)[7] * (f_7)[i] +
                                        (distribution_x)[8] * (f_8)[i]) / (h)[i];

                                (v)[i] = ((distribution_y)[0] * (f_0)[i] +
                                        (distribution_y)[1] * (f_1)[i] +
                                        (distribution_y)[2] * (f_2)[i] +
                                        (distribution_y)[3] * (f_3)[i] +
                                        (distribution_y)[4] * (f_4)[i] +
                                        (distribution_y)[5] * (f_5)[i] +
                                        (distribution_y)[6] * (f_6)[i] +
                                        (distribution_y)[7] * (f_7)[i] +
                                        (distribution_y)[8] * (f_8)[i]) / (h)[i];
                            }
                            else
                            {
                                //TODO: better heuristics for reset of h: if negative -> 0, if positive but too small -> epsilon
                                (h)[i] = 0;

                                (u)[i] = 0;
                                (v)[i] = 0;
                            }
                            (h)[i] = MinModLimiter<tags::CPU>::value((h)[i]);


                        }

                        info.limits->unlock(lm_read_only);

                        data.f_0->unlock(lm_read_and_write);
                        data.f_1->unlock(lm_read_and_write);
                        data.f_2->unlock(lm_read_and_write);
                        data.f_3->unlock(lm_read_and_write);
                        data.f_4->unlock(lm_read_and_write);
                        data.f_5->unlock(lm_read_and_write);
                        data.f_6->unlock(lm_read_and_write);
                        data.f_7->unlock(lm_read_and_write);
                        data.f_8->unlock(lm_read_and_write);

                        data.h->unlock(lm_write_only);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        data.u->unlock(lm_write_only);
                        data.v->unlock(lm_write_only);
                    }
        };

    template<>
        struct ExtractionGrid<tags::GPU::CUDA, lbm_modes::DRY>
        {
            public:
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon);
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data, double epsilon);
        };

    template<typename Tag_>
        struct ExtractionGrid<Tag_, lbm_modes::DRY>
        {
            public:
                template <typename DT_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ epsilon);
        };

    template<>
        struct ExtractionGrid<tags::Cell, lbm_modes::DRY>
        {
            public:
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon)
                {
                    ExtractionGrid<tags::CPU, lbm_modes::DRY>::value(info, data, epsilon);
                }
        };

    template<>
        struct ExtractionGrid<tags::CPU::SSE, lbm_modes::DRY>
        {
            public:
                template <typename DT1_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data, DT1_ epsilon);
        };

    template<>
        struct ExtractionGrid<tags::CPU::Itanium, lbm_modes::DRY>
        {
            public:
                template <typename DT1_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data, DT1_ epsilon)
                {
                    ExtractionGrid<tags::CPU, lbm_modes::DRY>::value(info, data, epsilon);
                }
        };

    template<>
        struct ExtractionGrid<tags::CPU, lbm_modes::WET>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, HONEI_UNUSED DT_ epsilon)
                    {
                        CONTEXT("When extracting h, u and v:");

                        //set f to t_temp
                        DenseVector<DT_> * swap;
                        swap = data.f_0;
                        data.f_0 = data.f_temp_0;
                        data.f_temp_0 = swap;
                        swap = data.f_1;
                        data.f_1 = data.f_temp_1;
                        data.f_temp_1 = swap;
                        swap = data.f_2;
                        data.f_2 = data.f_temp_2;
                        data.f_temp_2 = swap;
                        swap = data.f_3;
                        data.f_3 = data.f_temp_3;
                        data.f_temp_3 = swap;
                        swap = data.f_4;
                        data.f_4 = data.f_temp_4;
                        data.f_temp_4 = swap;
                        swap = data.f_5;
                        data.f_5 = data.f_temp_5;
                        data.f_temp_5 = swap;
                        swap = data.f_6;
                        data.f_6 = data.f_temp_6;
                        data.f_temp_6 = swap;
                        swap = data.f_7;
                        data.f_7 = data.f_temp_7;
                        data.f_temp_7 = swap;
                        swap = data.f_8;
                        data.f_8 = data.f_temp_8;
                        data.f_temp_8 = swap;

                        info.limits->lock(lm_read_only);

                        data.f_0->lock(lm_read_and_write);
                        data.f_1->lock(lm_read_and_write);
                        data.f_2->lock(lm_read_and_write);
                        data.f_3->lock(lm_read_and_write);
                        data.f_4->lock(lm_read_and_write);
                        data.f_5->lock(lm_read_and_write);
                        data.f_6->lock(lm_read_and_write);
                        data.f_7->lock(lm_read_and_write);
                        data.f_8->lock(lm_read_and_write);

                        data.h->lock(lm_write_only);

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        data.u->lock(lm_write_only);
                        data.v->lock(lm_write_only);

                        for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                        {

                            //accumulate
                            (*data.h)[i] = (*data.f_0)[i] +
                                (*data.f_1)[i] +
                                (*data.f_2)[i] +
                                (*data.f_3)[i] +
                                (*data.f_4)[i] +
                                (*data.f_5)[i] +
                                (*data.f_6)[i] +
                                (*data.f_7)[i] +
                                (*data.f_8)[i];

                            (*data.u)[i] = ((*data.distribution_x)[0] * (*data.f_0)[i] +
                                    (*data.distribution_x)[1] * (*data.f_1)[i] +
                                    (*data.distribution_x)[2] * (*data.f_2)[i] +
                                    (*data.distribution_x)[3] * (*data.f_3)[i] +
                                    (*data.distribution_x)[4] * (*data.f_4)[i] +
                                    (*data.distribution_x)[5] * (*data.f_5)[i] +
                                    (*data.distribution_x)[6] * (*data.f_6)[i] +
                                    (*data.distribution_x)[7] * (*data.f_7)[i] +
                                    (*data.distribution_x)[8] * (*data.f_8)[i]) / (*data.h)[i];

                            (*data.v)[i] = ((*data.distribution_y)[0] * (*data.f_0)[i] +
                                    (*data.distribution_y)[1] * (*data.f_1)[i] +
                                    (*data.distribution_y)[2] * (*data.f_2)[i] +
                                    (*data.distribution_y)[3] * (*data.f_3)[i] +
                                    (*data.distribution_y)[4] * (*data.f_4)[i] +
                                    (*data.distribution_y)[5] * (*data.f_5)[i] +
                                    (*data.distribution_y)[6] * (*data.f_6)[i] +
                                    (*data.distribution_y)[7] * (*data.f_7)[i] +
                                    (*data.distribution_y)[8] * (*data.f_8)[i]) / (*data.h)[i];


                        }

                        info.limits->unlock(lm_read_only);

                        data.f_0->unlock(lm_read_and_write);
                        data.f_1->unlock(lm_read_and_write);
                        data.f_2->unlock(lm_read_and_write);
                        data.f_3->unlock(lm_read_and_write);
                        data.f_4->unlock(lm_read_and_write);
                        data.f_5->unlock(lm_read_and_write);
                        data.f_6->unlock(lm_read_and_write);
                        data.f_7->unlock(lm_read_and_write);
                        data.f_8->unlock(lm_read_and_write);

                        data.h->unlock(lm_write_only);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        data.u->unlock(lm_write_only);
                        data.v->unlock(lm_write_only);
                    }

                template<typename DT1_>
                    static inline BenchmarkInfo get_benchmark_info(HONEI_UNUSED PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                    {
                        BenchmarkInfo result;
                        result.flops = data->h->size() * 44;
                        result.load = data->h->size() * (5 * 9 + 2) * sizeof(DT1_);
                        result.store = data->h->size() * 3 * sizeof(DT1_);
                        result.size.push_back(data->h->size());
                        return result;
                    }
        };

    template<>
        struct ExtractionGrid<tags::CPU::Generic, lbm_modes::WET>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, HONEI_UNUSED DT_ epsilon)
                    {
                        CONTEXT("When extracting h, u and v:");

                        //set f to t_temp
                        DenseVector<DT_> * swap;
                        swap = data.f_0;
                        data.f_0 = data.f_temp_0;
                        data.f_temp_0 = swap;
                        swap = data.f_1;
                        data.f_1 = data.f_temp_1;
                        data.f_temp_1 = swap;
                        swap = data.f_2;
                        data.f_2 = data.f_temp_2;
                        data.f_temp_2 = swap;
                        swap = data.f_3;
                        data.f_3 = data.f_temp_3;
                        data.f_temp_3 = swap;
                        swap = data.f_4;
                        data.f_4 = data.f_temp_4;
                        data.f_temp_4 = swap;
                        swap = data.f_5;
                        data.f_5 = data.f_temp_5;
                        data.f_temp_5 = swap;
                        swap = data.f_6;
                        data.f_6 = data.f_temp_6;
                        data.f_temp_6 = swap;
                        swap = data.f_7;
                        data.f_7 = data.f_temp_7;
                        data.f_temp_7 = swap;
                        swap = data.f_8;
                        data.f_8 = data.f_temp_8;
                        data.f_temp_8 = swap;

                        info.limits->lock(lm_read_only);

                        data.f_0->lock(lm_read_and_write);
                        data.f_1->lock(lm_read_and_write);
                        data.f_2->lock(lm_read_and_write);
                        data.f_3->lock(lm_read_and_write);
                        data.f_4->lock(lm_read_and_write);
                        data.f_5->lock(lm_read_and_write);
                        data.f_6->lock(lm_read_and_write);
                        data.f_7->lock(lm_read_and_write);
                        data.f_8->lock(lm_read_and_write);

                        data.h->lock(lm_write_only);

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        data.u->lock(lm_write_only);
                        data.v->lock(lm_write_only);

                        const unsigned long * const limits(info.limits->elements());
                        DT_ * f_0(data.f_0->elements());
                        DT_ * f_1(data.f_1->elements());
                        DT_ * f_2(data.f_2->elements());
                        DT_ * f_3(data.f_3->elements());
                        DT_ * f_4(data.f_4->elements());
                        DT_ * f_5(data.f_5->elements());
                        DT_ * f_6(data.f_6->elements());
                        DT_ * f_7(data.f_7->elements());
                        DT_ * f_8(data.f_8->elements());

                        DT_ * h(data.h->elements());
                        DT_ * u(data.u->elements());
                        DT_ * v(data.v->elements());
                        const DT_ * const distribution_x(data.distribution_x->elements());
                        const DT_ * const distribution_y(data.distribution_y->elements());

                        const unsigned long start(limits[0]);
                        const unsigned long end(limits[info.limits->size() - 1]);
                        for(unsigned long i(start); i < end; ++i)
                        {

                            //accumulate
                            (h)[i] = (f_0)[i] +
                                (f_1)[i] +
                                (f_2)[i] +
                                (f_3)[i] +
                                (f_4)[i] +
                                (f_5)[i] +
                                (f_6)[i] +
                                (f_7)[i] +
                                (f_8)[i];
                        }

                        for(unsigned long i(start); i < end; ++i)
                        {
                            (u)[i] = ((distribution_x)[0] * (f_0)[i] +
                                    (distribution_x)[1] * (f_1)[i] +
                                    (distribution_x)[2] * (f_2)[i] +
                                    (distribution_x)[3] * (f_3)[i] +
                                    (distribution_x)[4] * (f_4)[i] +
                                    (distribution_x)[5] * (f_5)[i] +
                                    (distribution_x)[6] * (f_6)[i] +
                                    (distribution_x)[7] * (f_7)[i] +
                                    (distribution_x)[8] * (f_8)[i]) / (h)[i];
                        }

                        for(unsigned long i(start); i < end; ++i)
                        {
                            (v)[i] = ((distribution_y)[0] * (f_0)[i] +
                                    (distribution_y)[1] * (f_1)[i] +
                                    (distribution_y)[2] * (f_2)[i] +
                                    (distribution_y)[3] * (f_3)[i] +
                                    (distribution_y)[4] * (f_4)[i] +
                                    (distribution_y)[5] * (f_5)[i] +
                                    (distribution_y)[6] * (f_6)[i] +
                                    (distribution_y)[7] * (f_7)[i] +
                                    (distribution_y)[8] * (f_8)[i]) / (h)[i];
                        }

                        info.limits->unlock(lm_read_only);

                        data.f_0->unlock(lm_read_and_write);
                        data.f_1->unlock(lm_read_and_write);
                        data.f_2->unlock(lm_read_and_write);
                        data.f_3->unlock(lm_read_and_write);
                        data.f_4->unlock(lm_read_and_write);
                        data.f_5->unlock(lm_read_and_write);
                        data.f_6->unlock(lm_read_and_write);
                        data.f_7->unlock(lm_read_and_write);
                        data.f_8->unlock(lm_read_and_write);

                        data.h->unlock(lm_write_only);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        data.u->unlock(lm_write_only);
                        data.v->unlock(lm_write_only);
                    }
        };

    template<>
        struct ExtractionGrid<tags::GPU::CUDA, lbm_modes::WET>
        {
            public:
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon);
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data, double epsilon);
        };

    template<typename Tag_>
        struct ExtractionGrid<Tag_, lbm_modes::WET>
        {
            public:
                template <typename DT_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ epsilon);
        };

    template<>
        struct ExtractionGrid<tags::CPU::SSE, lbm_modes::WET>
        {
            public:
                template <typename DT1_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data, DT1_ epsilon);
        };

    template<>
        struct ExtractionGrid<tags::Cell, lbm_modes::WET>
        {
            public:
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon);
        };

    template<>
        struct ExtractionGrid<tags::CPU::Itanium, lbm_modes::WET>
        {
            public:
                template <typename DT1_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data, DT1_ epsilon);
        };
}
#endif
