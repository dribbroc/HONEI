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

namespace quantities
{
    class HEIGHT;
    class VELOCITY_X;
    class VELOCITY_Y;
    class DENSITY; //used for NAVSTO later
}

namespace honei
{
    template<typename Tag_, typename App_, typename Quantity_>
        class Extraction
        {
        };

    template<>
        class Extraction<tags::CPU, lbm_applications::LABSWE , quantities::HEIGHT>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
                    {
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
                        }
                    }
        };


    template<>
        class Extraction<tags::CPU, lbm_applications::LABSWE , quantities::VELOCITY_X>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
                    {
                        for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                        {
                            //accumulate
                            (*data.u)[i] = ((*data.distribution_x)[0] * (*data.f_0)[i] +
                                    (*data.distribution_x)[1] * (*data.f_1)[i] +
                                    (*data.distribution_x)[2] * (*data.f_2)[i] +
                                    (*data.distribution_x)[3] * (*data.f_3)[i] +
                                    (*data.distribution_x)[4] * (*data.f_4)[i] +
                                    (*data.distribution_x)[5] * (*data.f_5)[i] +
                                    (*data.distribution_x)[6] * (*data.f_6)[i] +
                                    (*data.distribution_x)[7] * (*data.f_7)[i] +
                                    (*data.distribution_x)[8] * (*data.f_8)[i]) / (*data.h)[i];
                        }
                    }
        };


    template<>
        class Extraction<tags::CPU, lbm_applications::LABSWE , quantities::VELOCITY_Y>
        {
            public:
                template<typename DT_>
                    static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
                    {
                        for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                        {
                            //accumulate
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
                    }
        };
}
#endif
