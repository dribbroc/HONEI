/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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


using namespace honei;
using namespace lbm;

namespace honei
{
    template<typename Tag_, typename App_>
        struct Extraction
        {
        };

    template<>
        struct Extraction<tags::CPU, lbm_applications::LABSWE>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
                    {
                        CONTEXT("When extracting h, u and v:");

                        info.limits->lock(lm_read_only);

                        data.f_temp_0->lock(lm_read_only);
                        data.f_temp_1->lock(lm_read_only);
                        data.f_temp_2->lock(lm_read_only);
                        data.f_temp_3->lock(lm_read_only);
                        data.f_temp_4->lock(lm_read_only);
                        data.f_temp_5->lock(lm_read_only);
                        data.f_temp_6->lock(lm_read_only);
                        data.f_temp_7->lock(lm_read_only);
                        data.f_temp_8->lock(lm_read_only);

                        data.f_0->lock(lm_read_and_write);
                        data.f_1->lock(lm_read_and_write);
                        data.f_2->lock(lm_read_and_write);
                        data.f_3->lock(lm_read_and_write);
                        data.f_4->lock(lm_read_and_write);
                        data.f_5->lock(lm_read_and_write);
                        data.f_6->lock(lm_read_and_write);
                        data.f_7->lock(lm_read_and_write);
                        data.f_8->lock(lm_read_and_write);

                        data.h->lock(lm_write_only); // in this case: write before read

                        data.distribution_x->lock(lm_read_only);
                        data.distribution_y->lock(lm_read_only);

                        data.u->lock(lm_write_only);
                        data.v->lock(lm_write_only);

                        for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                        {
                            //set f to t_temp
                            (*data.f_0)[i] = (*data.f_temp_0)[i];
                            (*data.f_1)[i] = (*data.f_temp_1)[i];
                            (*data.f_2)[i] = (*data.f_temp_2)[i];
                            (*data.f_3)[i] = (*data.f_temp_3)[i];
                            (*data.f_4)[i] = (*data.f_temp_4)[i];
                            (*data.f_5)[i] = (*data.f_temp_5)[i];
                            (*data.f_6)[i] = (*data.f_temp_6)[i];
                            (*data.f_7)[i] = (*data.f_temp_7)[i];
                            (*data.f_8)[i] = (*data.f_temp_8)[i];

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

                        data.f_temp_0->unlock(lm_read_only);
                        data.f_temp_1->unlock(lm_read_only);
                        data.f_temp_2->unlock(lm_read_only);
                        data.f_temp_3->unlock(lm_read_only);
                        data.f_temp_4->unlock(lm_read_only);
                        data.f_temp_5->unlock(lm_read_only);
                        data.f_temp_6->unlock(lm_read_only);
                        data.f_temp_7->unlock(lm_read_only);
                        data.f_temp_8->unlock(lm_read_only);

                        data.f_0->unlock(lm_read_and_write);
                        data.f_1->unlock(lm_read_and_write);
                        data.f_2->unlock(lm_read_and_write);
                        data.f_3->unlock(lm_read_and_write);
                        data.f_4->unlock(lm_read_and_write);
                        data.f_5->unlock(lm_read_and_write);
                        data.f_6->unlock(lm_read_and_write);
                        data.f_7->unlock(lm_read_and_write);
                        data.f_8->unlock(lm_read_and_write);

                        data.h->unlock(lm_read_and_write);

                        data.distribution_x->unlock(lm_read_only);
                        data.distribution_y->unlock(lm_read_only);

                        data.u->unlock(lm_write_only);
                        data.v->unlock(lm_write_only);
                    }
        };

    template<>
        struct Extraction<tags::GPU::CUDA, lbm_applications::LABSWE>
        {
            public:
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data);
        };
}
#endif
