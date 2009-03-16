/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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


#ifndef LBM_GUARD_FORCE_GRID_HH
#define LBM_GUARD_FORCE_GRID_HH 1

/**
 * \file
 * Implementation of source modules used by  LBM - (SWE) solvers using a
 * PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/lbm/grid.hh>
#include <honei/la/dense_vector.hh>
#include <honei/util/benchmark_info.hh>
#include <cmath>
using namespace honei::lbm;

namespace honei
{
   template <typename Tag_,
              typename App_,
              typename SourceType_,
              typename SourceScheme_>
    struct ForceGrid
    {
    };

    /**
     * \brief Centred scheme bed slope force term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <>
    struct ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>
    {
        /**
         * \name Force term.
         *
         * \brief Computes a simple source term value.
         *
         * \param data Our packed grid data object.
         * \param info Our info object.
         * \param d_x Our delta x.
         * \param d_y Our delta y.
         * \param d_t Our delta t.
         * \param g The gravitational constant to be used.
         *
         */
        template<typename DT1_, typename DT2_>
            static void value(PackedGridData<D2Q9, DT1_> & data, PackedGridInfo<D2Q9> & info, DT2_ g, DT2_ d_x, DT2_ d_y, DT2_ d_t, DT2_ manning)
            {
                CONTEXT("When computing LABSWE force term:");

                info.dir_index_1->lock(lm_read_only);
                info.dir_index_2->lock(lm_read_only);
                info.dir_index_3->lock(lm_read_only);
                info.dir_index_4->lock(lm_read_only);
                info.dir_index_5->lock(lm_read_only);
                info.dir_index_6->lock(lm_read_only);
                info.dir_index_7->lock(lm_read_only);
                info.dir_index_8->lock(lm_read_only);

                info.dir_1->lock(lm_read_only);
                info.dir_2->lock(lm_read_only);
                info.dir_3->lock(lm_read_only);
                info.dir_4->lock(lm_read_only);
                info.dir_5->lock(lm_read_only);
                info.dir_6->lock(lm_read_only);
                info.dir_7->lock(lm_read_only);
                info.dir_8->lock(lm_read_only);

                data.h->lock(lm_read_only);
                data.b->lock(lm_read_only);
                data.distribution_x->lock(lm_read_only);
                data.distribution_y->lock(lm_read_only);

                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

                //precompute constants
                DT1_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
                DT1_ gravity_multiplier(-g);

                //-----------alpha = 1 ----------------------------------------------------------------------------------------------
                unsigned long size_1(data.f_temp_1->size());
                DenseVector<DT1_> temp_1(size_1, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_1[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[1];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_1[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_1)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_1[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp_1 to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_1)[i] += temp_1[i];
                }
                //repeat for y direction
                //NOTHING TO BE DONE HERE


                //-----------alpha = 2 ----------------------------------------------------------------------------------------------

                unsigned long size_2(data.f_temp_2->size());
                DenseVector<DT1_> temp_2(size_2, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[2];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_2)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x) (still use direction 1 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp_2 to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_2)[i] += temp_2[i];
                }

                //REPEAT FOR Y

                DenseVector<DT1_> temp_2_y(size_2, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[2];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_2)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_2_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp_2_y to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_2)[i] += temp_2_y[i];
                }


                //-----------alpha = 3 ----------------------------------------------------------------------------------------------

                // Y DIRECTION ONLY

                unsigned long size_3(data.f_temp_3->size());
                DenseVector<DT1_> temp_3_y(size_3, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_3_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[3];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_3_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_3)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_3_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp_3_y to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_3)[i] += temp_3_y[i];
                }

                //-----------alpha = 4 ----------------------------------------------------------------------------------------------

                unsigned long size_4(data.f_temp_4->size());
                DenseVector<DT1_> temp_4(size_4, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[4];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_4)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x) (still use direction 1 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_4)[i] += temp_4[i];
                }

                //REPEAT FOR Y

                DenseVector<DT1_> temp_4_y(size_4, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[4];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_4)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_4_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_4)[i] += temp_4_y[i];
                }

                //-----------alpha = 5 ----------------------------------------------------------------------------------------------
                //X ONLY
                unsigned long size_5(data.f_temp_5->size());
                DenseVector<DT1_> temp_5(size_5, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                    {
                        temp_5[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[5];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                    {
                        temp_5[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_5)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x) (still use direction 1 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_5[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_5)[i] += temp_5[i];
                }

                //-----------alpha = 6 ----------------------------------------------------------------------------------------------

                unsigned long size_6(data.f_temp_6->size());
                DenseVector<DT1_> temp_6(size_6, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[6];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_6)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x) (still use direction 1 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_6)[i] += temp_6[i];
                }

                //REPEAT FOR Y

                DenseVector<DT1_> temp_6_y(size_6, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[6];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_6)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_6_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_6)[i] += temp_6_y[i];
                }

                //-----------alpha = 7 ----------------------------------------------------------------------------------------------

                //Y ONLY

                unsigned long size_7(data.f_temp_7->size());
                DenseVector<DT1_> temp_7_y(size_7, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                    {
                        temp_7_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[7];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                    {
                        temp_7_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_7)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_7_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_7)[i] += temp_7_y[i];
                }

                //-----------alpha = 8 ----------------------------------------------------------------------------------------------

                unsigned long size_8(data.f_temp_8->size());
                DenseVector<DT1_> temp_8(size_8, DT1_(0));

                //set up temp (x)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8[i] = force_multiplier * gravity_multiplier * (*data.distribution_x)[8];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_8)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (x) (still use direction 1 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8[i] *= ((*data.b)[(*info.dir_1)[half] + offset] - (*data.b)[i]) / (d_x);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_8)[i] += temp_8[i];
                }

                //REPEAT FOR Y

                DenseVector<DT1_> temp_8_y(size_8, DT1_(0));

                //set up temp (y)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8_y[i] = force_multiplier * gravity_multiplier * (*data.distribution_y)[8];
                    }
                }
                //multiply temp by interpolation
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8_y[i] *= (((*data.h)[i]) + ((*data.h)[(*info.dir_8)[half] + offset])) / (DT1_(2));
                    }
                }
                //multiply temp by derivative (y) (use direction 3 due to forward differences)
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        temp_8_y[i] *= ((*data.b)[(*info.dir_3)[half] + offset] - (*data.b)[i]) / (d_y);
                    }
                }
                //add temp to distribution
                for(unsigned long i((*info.limits)[0]) ; i < (*info.limits)[info.limits->size() -1] ; ++i)
                {
                    (*data.f_temp_8)[i] += temp_8_y[i];
                }

                info.dir_index_1->unlock(lm_read_only);
                info.dir_index_2->unlock(lm_read_only);
                info.dir_index_3->unlock(lm_read_only);
                info.dir_index_4->unlock(lm_read_only);
                info.dir_index_5->unlock(lm_read_only);
                info.dir_index_6->unlock(lm_read_only);
                info.dir_index_7->unlock(lm_read_only);
                info.dir_index_8->unlock(lm_read_only);

                info.dir_1->unlock(lm_read_only);
                info.dir_2->unlock(lm_read_only);
                info.dir_3->unlock(lm_read_only);
                info.dir_4->unlock(lm_read_only);
                info.dir_5->unlock(lm_read_only);
                info.dir_6->unlock(lm_read_only);
                info.dir_7->unlock(lm_read_only);
                info.dir_8->unlock(lm_read_only);

                data.h->unlock(lm_read_only);
                data.b->unlock(lm_read_only);
                data.distribution_x->unlock(lm_read_only);
                data.distribution_y->unlock(lm_read_only);

                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }

            template<typename DT1_>
                static inline BenchmarkInfo get_benchmark_info(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                {
                    BenchmarkInfo result;
                    result.flops = data->h->size() * 60;
                    result.load = data->h->size() * 40 * sizeof(DT1_) + 3 * sizeof(DT1_);
                    result.store = data->h->size() * 9 * sizeof(DT1_);
                    result.size.push_back(data->h->size());
                    return result;
                }
    };


    /**
     * \brief Centred scheme bed friction force term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <>
    struct ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>
    {
        /**
         * \name Force term.
         *
         * \brief Computes a simple source term value.
         *
         * \param data Our packed grid data object.
         * \param info Our info object.
         * \param d_x Our delta x.
         * \param d_y Our delta y.
         * \param d_t Our delta t.
         * \param g The gravitational constant to be used.
         *
         */
        template<typename DT1_, typename DT2_>
            static void value(PackedGridData<D2Q9, DT1_> & data, PackedGridInfo<D2Q9> & info, DT2_ g, DT2_ d_x, DT2_ d_y, DT2_ d_t, DT2_ manning_const)
            {
                CONTEXT("When computing LABSWE force term:");

                info.dir_index_1->lock(lm_read_only);
                info.dir_index_2->lock(lm_read_only);
                info.dir_index_3->lock(lm_read_only);
                info.dir_index_4->lock(lm_read_only);
                info.dir_index_5->lock(lm_read_only);
                info.dir_index_6->lock(lm_read_only);
                info.dir_index_7->lock(lm_read_only);
                info.dir_index_8->lock(lm_read_only);

                info.dir_1->lock(lm_read_only);
                info.dir_2->lock(lm_read_only);
                info.dir_3->lock(lm_read_only);
                info.dir_4->lock(lm_read_only);
                info.dir_5->lock(lm_read_only);
                info.dir_6->lock(lm_read_only);
                info.dir_7->lock(lm_read_only);
                info.dir_8->lock(lm_read_only);

                data.h->lock(lm_read_only);
                data.b->lock(lm_read_only);
                data.distribution_x->lock(lm_read_only);
                data.distribution_y->lock(lm_read_only);

                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

                //precompute constants
                DT1_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));

                //-----------alpha = 1 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()) )
                        {
                            (*data.f_temp_1)[i] -= force_multiplier * g * (*data.distribution_x)[1] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }


                //-----------alpha = 2 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_2)[i] -= force_multiplier * g * (*data.distribution_x)[2] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_2)[i] -= force_multiplier * g * (*data.distribution_y)[2] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }


                //-----------alpha = 3 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_3)[i] -= force_multiplier * g * (*data.distribution_y)[3] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                //-----------alpha = 4 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_4)[i] -= force_multiplier * g * (*data.distribution_x)[4] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_4)[i] -= force_multiplier * g * (*data.distribution_y)[4] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }


                //-----------alpha = 5 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_5)[i] -= force_multiplier * g * (*data.distribution_x)[5] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                //-----------alpha = 6 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_6)[i] -= force_multiplier * g * (*data.distribution_x)[6] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_6)[i] -= force_multiplier * g * (*data.distribution_y)[6] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                //-----------alpha = 7 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_7)[i] -= force_multiplier * g * (*data.distribution_y)[7] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }


                //-----------alpha = 8 ----------------------------------------------------------------------------------------------

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_8)[i] -= force_multiplier * g * (*data.distribution_x)[8] * (manning_const * manning_const) *
                                (*data.u)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        if ( (pow((*data.h)[i], DT2_(1./3.))) > std::numeric_limits<DT2_>::epsilon() || (pow((*data.h)[i], DT2_(1./3.))) < DT2_(-std::numeric_limits<DT2_>::epsilon()))
                        {
                            (*data.f_temp_8)[i] -= force_multiplier * g * (*data.distribution_y)[8] * (manning_const * manning_const) *
                                (*data.v)[i] * sqrt((*data.u)[i] * (*data.u)[i] + (*data.v)[i] * (*data.v)[i]) / (pow((*data.h)[i], DT2_(1./3.)));
                        }
                    }
                }

                info.dir_index_1->unlock(lm_read_only);
                info.dir_index_2->unlock(lm_read_only);
                info.dir_index_3->unlock(lm_read_only);
                info.dir_index_4->unlock(lm_read_only);
                info.dir_index_5->unlock(lm_read_only);
                info.dir_index_6->unlock(lm_read_only);
                info.dir_index_7->unlock(lm_read_only);
                info.dir_index_8->unlock(lm_read_only);

                info.dir_1->unlock(lm_read_only);
                info.dir_2->unlock(lm_read_only);
                info.dir_3->unlock(lm_read_only);
                info.dir_4->unlock(lm_read_only);
                info.dir_5->unlock(lm_read_only);
                info.dir_6->unlock(lm_read_only);
                info.dir_7->unlock(lm_read_only);
                info.dir_8->unlock(lm_read_only);

                data.h->unlock(lm_read_only);
                data.b->unlock(lm_read_only);
                data.distribution_x->unlock(lm_read_only);
                data.distribution_y->unlock(lm_read_only);

                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }

            template<typename DT1_>
            static inline BenchmarkInfo get_benchmark_info(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                {
                    // \todo Insert real data
                    BenchmarkInfo result;
                    result.flops = data->h->size() * 44 + data->h->size() * 88;
                    result.load = data->h->size() * 40 * sizeof(DT1_) + 5 * sizeof(DT1_);
                    result.store = data->h->size() * 9 * sizeof(DT1_);
                    result.size.push_back(data->h->size());
                    return result;
                }
    };

    template <typename Tag_>
    struct ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FULL>
    {
        template<typename DT1_, typename DT2_>
            static void value(PackedGridData<D2Q9, DT1_> & data, PackedGridInfo<D2Q9> & info, DT2_ g, DT2_ d_x, DT2_ d_y, DT2_ d_t, DT2_ manning)
            {
                ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(data, info, g, d_x, d_y, d_t, manning);
                ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(data, info, g, d_x, d_y, d_t, manning);
            }

            template<typename DT1_>
                static inline BenchmarkInfo get_benchmark_info(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                {
                    BenchmarkInfo result;
                    BenchmarkInfo f1(ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::get_benchmark_info(info, data));
                    result + f1;
                    BenchmarkInfo f2(ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::get_benchmark_info(info, data));
                    result + f2;
                    result.size.push_back(data->h->size());
                    return result;
                }
    };

    template <>
        struct ForceGrid<tags::GPU::CUDA, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>
        {
            static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t, float manning)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(data, info, g, d_x, d_y, d_t, manning);
            }
        };

    template <>
        struct ForceGrid<tags::GPU::CUDA, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>
        {
            static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t, float m)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(data, info, g, d_x, d_y, d_t, m);
            }
        };

    template <>
        struct ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>
        {
            static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t, float manning)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(data, info, g, d_x, d_y, d_t, manning);
            }
            static void value(PackedGridData<D2Q9, double> & data, PackedGridInfo<D2Q9> & info, double g, double d_x, double d_y, double d_t, double manning)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(data, info, g, d_x, d_y, d_t, manning);
            }
        };

    template <>
        struct ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>
        {
            static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t, float m)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(data, info, g, d_x, d_y, d_t, m);
            }
            static void value(PackedGridData<D2Q9, double> & data, PackedGridInfo<D2Q9> & info, double g, double d_x, double d_y, double d_t, double m)
            {
                ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(data, info, g, d_x, d_y, d_t, m);
            }
        };
}
#endif
