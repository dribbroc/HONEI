/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <apryde@gmx.de>
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

#ifndef LBM_GUARD_EXTRACTION_GRID_POLLUTANT_HH
#define LBM_GUARD_EXTRACTION_GRID_POLLUTANT_HH 1


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
        struct ExtractionGridPollutant
        {
        };

    template<>
        struct ExtractionGridPollutant<tags::CPU::Generic, lbm_modes::DRY>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data_flow, PackedGridData<D2Q9, DT_> & data_poll, HONEI_UNUSED DT_ epsilon)
                    {
                        CONTEXT("When extracting c only:");

                        //set f to t_temp
                        DenseVector<DT_> * swap;
                        swap = data_poll.f_0;
                        data_poll.f_0 = data_poll.f_temp_0;
                        data_poll.f_temp_0 = swap;
                        swap = data_poll.f_1;
                        data_poll.f_1 = data_poll.f_temp_1;
                        data_poll.f_temp_1 = swap;
                        swap = data_poll.f_2;
                        data_poll.f_2 = data_poll.f_temp_2;
                        data_poll.f_temp_2 = swap;
                        swap = data_poll.f_3;
                        data_poll.f_3 = data_poll.f_temp_3;
                        data_poll.f_temp_3 = swap;
                        swap = data_poll.f_4;
                        data_poll.f_4 = data_poll.f_temp_4;
                        data_poll.f_temp_4 = swap;
                        swap = data_poll.f_5;
                        data_poll.f_5 = data_poll.f_temp_5;
                        data_poll.f_temp_5 = swap;
                        swap = data_poll.f_6;
                        data_poll.f_6 = data_poll.f_temp_6;
                        data_poll.f_temp_6 = swap;
                        swap = data_poll.f_7;
                        data_poll.f_7 = data_poll.f_temp_7;
                        data_poll.f_temp_7 = swap;
                        swap = data_poll.f_8;
                        data_poll.f_8 = data_poll.f_temp_8;
                        data_poll.f_temp_8 = swap;

                        info.limits->lock(lm_read_only);

                        data_poll.f_0->lock(lm_read_and_write);
                        data_poll.f_1->lock(lm_read_and_write);
                        data_poll.f_2->lock(lm_read_and_write);
                        data_poll.f_3->lock(lm_read_and_write);
                        data_poll.f_4->lock(lm_read_and_write);
                        data_poll.f_5->lock(lm_read_and_write);
                        data_poll.f_6->lock(lm_read_and_write);
                        data_poll.f_7->lock(lm_read_and_write);
                        data_poll.f_8->lock(lm_read_and_write);

                        data_flow.h->lock(lm_read_only);
                        data_poll.h->lock(lm_write_only);

                        const unsigned long * const limits(info.limits->elements());
                        DT_ * f_0(data_poll.f_0->elements());
                        DT_ * f_1(data_poll.f_1->elements());
                        DT_ * f_2(data_poll.f_2->elements());
                        DT_ * f_3(data_poll.f_3->elements());
                        DT_ * f_4(data_poll.f_4->elements());
                        DT_ * f_5(data_poll.f_5->elements());
                        DT_ * f_6(data_poll.f_6->elements());
                        DT_ * f_7(data_poll.f_7->elements());
                        DT_ * f_8(data_poll.f_8->elements());

                        DT_ * h(data_flow.h->elements());
                        DT_ * c(data_poll.h->elements());

                        //DT_ lax_upper(epsilon /*10e-5 std::numeric_limits<DT_>::epsilon() * 10e7*/);
                        //DT_ lax_lower(-lax_upper);

                        const unsigned long start(limits[0]);
                        const unsigned long end(limits[info.limits->size() - 1]);
                        for(unsigned long i(start); i < end; ++i)
                        {

                            //accumulate
                            (c)[i] =( (f_0)[i] +
                                (f_1)[i] +
                                (f_2)[i] +
                                (f_3)[i] +
                                (f_4)[i] +
                                (f_5)[i] +
                                (f_6)[i] +
                                (f_7)[i] +
                                (f_8)[i] ) / h[i];

                            //if((h)[i] < lax_lower || (h)[i] > lax_upper)
                            //{
                            //}
                            //else
                            //{
                                //TODO: better heuristics for reset of h: if negative -> 0, if positive but too small -> epsilon
                                //(c)[i] = 0;
                            //}
                            //(h)[i] = MinModLimiter<tags::CPU>::value((h)[i]);

                        }

                        info.limits->unlock(lm_read_only);

                        data_poll.f_0->unlock(lm_read_and_write);
                        data_poll.f_1->unlock(lm_read_and_write);
                        data_poll.f_2->unlock(lm_read_and_write);
                        data_poll.f_3->unlock(lm_read_and_write);
                        data_poll.f_4->unlock(lm_read_and_write);
                        data_poll.f_5->unlock(lm_read_and_write);
                        data_poll.f_6->unlock(lm_read_and_write);
                        data_poll.f_7->unlock(lm_read_and_write);
                        data_poll.f_8->unlock(lm_read_and_write);

                        data_flow.h->unlock(lm_read_only);
                        data_poll.h->unlock(lm_write_only);

                    }
        };



}


#endif
