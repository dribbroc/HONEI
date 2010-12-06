/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/lbm/force_grid.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;



template <typename DT_>
void init_temp(const unsigned long * dir, const unsigned long * dir_index, unsigned long dir_index_size,
        DT_ * temp, const DT_ * h, DT_ dist, DT_ ftg)
{
    DT_ fac(ftg * dist);
    for (unsigned long begin(0), half(0) ; begin < dir_index_size - 1; begin+=2, ++half)
    {
        for (unsigned long i(dir_index[begin]), offset(0) ; i < dir_index[begin + 1] ; ++i, ++offset)
        {
            temp[i] = fac * ((h[i] + (h[dir[half] + offset])) / (DT_(2)));
        }
    }
}

template <typename DT_>
void derive_temp(const unsigned long * dir, const unsigned long * dir_index, unsigned long dir_index_size,
        DT_ * temp, const DT_ * b, DT_ dxy)
{
    for (unsigned long begin(0), half(0) ; begin < dir_index_size - 1; begin+=2, ++half)
    {
        for (unsigned long i(dir_index[begin]), offset(0) ; i < dir_index[begin + 1] ; ++i, ++offset)
        {
            temp[i] *= (b[dir[half] + offset] - b[i]) / (dxy);
        }
    }
}

template <typename DT_>
void add_temp(DT_ * f, const DT_ * temp, unsigned long start, unsigned long end)
{
    for (unsigned long i(start) ; i < end ; ++i)
    {
        f[i] += temp[i];
    }
}


template <typename DT_>
void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(
            PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ g, DT_ d_x, DT_ d_y, DT_ d_t, HONEI_UNUSED DT_ manning)
{
    CONTEXT("When computing LABSWE force term: SLOPE SSE");

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
    DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
    DT_ gravity_multiplier(-g);
    DT_ force_times_gravity(force_multiplier * gravity_multiplier);

    DenseVector<DT_> temp(*data.temp);

    //-----------alpha = 1 ----------------------------------------------------------------------------------------------

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[1], force_times_gravity);

    //multiply temp by derivative (x)
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_1->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);
    //repeat for y direction
    //NOTHING TO BE DONE HERE


    //-----------alpha = 2 ----------------------------------------------------------------------------------------------

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_2->elements(), info.dir_index_2->elements(), info.dir_index_2->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[2], force_times_gravity);
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_2->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //REPEAT FOR Y

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_2->elements(), info.dir_index_2->elements(), info.dir_index_2->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[2], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_2->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);


    //-----------alpha = 3 ----------------------------------------------------------------------------------------------

    // Y DIRECTION ONLY

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[3], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_3->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //-----------alpha = 4 ----------------------------------------------------------------------------------------------

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_4->elements(), info.dir_index_4->elements(), info.dir_index_4->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[4], force_times_gravity);
    //multiply temp by derivative (x) (still use direction 1 due to forward differences)
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_4->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //REPEAT FOR Y

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_4->elements(), info.dir_index_4->elements(), info.dir_index_4->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[4], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_4->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //-----------alpha = 5 ----------------------------------------------------------------------------------------------
    //X ONLY

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_5->elements(), info.dir_index_5->elements(), info.dir_index_5->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[5], force_times_gravity);
    //multiply temp by derivative (x) (still use direction 1 due to forward differences)
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_5->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //-----------alpha = 6 ----------------------------------------------------------------------------------------------

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_6->elements(), info.dir_index_6->elements(), info.dir_index_6->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[6], force_times_gravity);
    //multiply temp by derivative (x) (still use direction 1 due to forward differences)
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_6->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //REPEAT FOR Y

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_6->elements(), info.dir_index_6->elements(), info.dir_index_6->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[6], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_6->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //-----------alpha = 7 ----------------------------------------------------------------------------------------------

    //Y ONLY

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_7->elements(), info.dir_index_7->elements(), info.dir_index_7->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[7], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_7->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //-----------alpha = 8 ----------------------------------------------------------------------------------------------

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_8->elements(), info.dir_index_8->elements(), info.dir_index_8->size(), temp.elements(), data.h->elements(), (*data.distribution_x)[8], force_times_gravity);
    //multiply temp by derivative (x) (still use direction 1 due to forward differences)
    derive_temp(info.dir_1->elements(), info.dir_index_1->elements(), info.dir_index_1->size(), temp.elements(), data.b->elements(), d_x);
    add_temp(data.f_temp_8->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

    //REPEAT FOR Y

    fill<tags::CPU::SSE>(temp);
    init_temp(info.dir_8->elements(), info.dir_index_8->elements(), info.dir_index_8->size(), temp.elements(), data.h->elements(), (*data.distribution_y)[8], force_times_gravity);
    //multiply temp by derivative (y) (use direction 3 due to forward differences)
    derive_temp(info.dir_3->elements(), info.dir_index_3->elements(), info.dir_index_3->size(), temp.elements(), data.b->elements(), d_y);
    add_temp(data.f_temp_8->elements(), temp.elements(), (*info.limits)[0], (*info.limits)[info.limits->size() - 1]);

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

template void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, float> &, float, float, float, float, float);

template void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, double> &, double, double, double, double, double);



void force_friction(const unsigned long start, const unsigned long end, double * f_temp, const double * h,
        const double * u, const double * v, double force_multiplier, double dist, double manning_const_sq)
{
    sse::force_friction(start, end, f_temp, h, u, v, force_multiplier, dist, manning_const_sq);
}

inline void force_friction(const unsigned long start, const unsigned long end, float * f_temp, const float * h,
        const float * u, const float * v, float force_multiplier, float dist, float manning_const_sq)
{
    sse::force_friction(start, end, f_temp, h, u, v, force_multiplier, dist, manning_const_sq);
}

template <typename DT_>
void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ g, DT_ d_x, HONEI_UNUSED DT_ d_y, DT_ d_t, DT_ manning_const)
{
    CONTEXT("When computing LABSWE force term: FRICTION SSE");

    info.limits->lock(lm_read_only);

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
    DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
    force_multiplier *= g;
    DT_ manning_const_sq(manning_const * manning_const);

    //-----------alpha = 1 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_1->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[1], manning_const_sq);


    //-----------alpha = 2 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_2->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[2], manning_const_sq);

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_2->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[2], manning_const_sq);


    //-----------alpha = 3 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_3->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[3], manning_const_sq);

    //-----------alpha = 4 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_4->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[4], manning_const_sq);

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_4->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[4], manning_const_sq);


    //-----------alpha = 5 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_5->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[5], manning_const_sq);

    //-----------alpha = 6 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_6->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[6], manning_const_sq);

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_6->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[6], manning_const_sq);

    //-----------alpha = 7 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_7->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[7], manning_const_sq);


    //-----------alpha = 8 ----------------------------------------------------------------------------------------------

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_8->elements(), data.h->elements(),
            data.u->elements(), data.v->elements(), force_multiplier, (*data.distribution_x)[8], manning_const_sq);

    force_friction((*info.limits)[0], (*info.limits)[info.limits->size() - 1], data.f_temp_8->elements(), data.h->elements(),
            data.v->elements(), data.u->elements(), force_multiplier, (*data.distribution_y)[8], manning_const_sq);



    info.limits->unlock(lm_read_only);

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


template void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, float> &, float, float, float, float, float);

template void ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, double> &, double, double, double, double, double);
