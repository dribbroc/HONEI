/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/backends/itanium/operations.hh>


using namespace honei;

template <typename DT1_>
void EquilibriumDistributionGrid<tags::CPU::Itanium, lbm_applications::LABSWE>::value(DT1_ g, DT1_ e,
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data)
{
    CONTEXT("When computing LABSWE local equilibrium distribution function (Itanium):");

    info.limits->lock(lm_read_only);

    data.u->lock(lm_read_only);
    data.v->lock(lm_read_only);
    data.h->lock(lm_read_only);

    data.distribution_x->lock(lm_read_only);
    data.distribution_y->lock(lm_read_only);

    data.f_eq_0->lock(lm_write_only);
    data.f_eq_1->lock(lm_write_only);
    data.f_eq_2->lock(lm_write_only);
    data.f_eq_3->lock(lm_write_only);
    data.f_eq_4->lock(lm_write_only);
    data.f_eq_5->lock(lm_write_only);
    data.f_eq_6->lock(lm_write_only);
    data.f_eq_7->lock(lm_write_only);
    data.f_eq_8->lock(lm_write_only);

    unsigned long begin((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);


    itanium::eq_dist_grid_dir_0(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.f_eq_0->elements());

    itanium::eq_dist_grid_dir_odd(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_1->elements(), 1);

    itanium::eq_dist_grid_dir_odd(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_3->elements(), 3);

    itanium::eq_dist_grid_dir_odd(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_5->elements(), 5);

    itanium::eq_dist_grid_dir_odd(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_7->elements(), 7);

    itanium::eq_dist_grid_dir_even(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_2->elements(), 2);

    itanium::eq_dist_grid_dir_even(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_4->elements(), 4);

    itanium::eq_dist_grid_dir_even(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_6->elements(), 6);

    itanium::eq_dist_grid_dir_even(begin, end, g, e,
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.f_eq_8->elements(), 8);

    info.limits->unlock(lm_read_only);

    data.u->unlock(lm_read_only);
    data.v->unlock(lm_read_only);
    data.h->unlock(lm_read_only);

    data.distribution_x->unlock(lm_read_only);
    data.distribution_y->unlock(lm_read_only);

    data.f_eq_0->unlock(lm_write_only);
    data.f_eq_1->unlock(lm_write_only);
    data.f_eq_2->unlock(lm_write_only);
    data.f_eq_3->unlock(lm_write_only);
    data.f_eq_4->unlock(lm_write_only);
    data.f_eq_5->unlock(lm_write_only);
    data.f_eq_6->unlock(lm_write_only);
    data.f_eq_7->unlock(lm_write_only);
    data.f_eq_8->unlock(lm_write_only);
}

template void EquilibriumDistributionGrid<tags::CPU::Itanium, lbm_applications::LABSWE>::value<float>(float, float,
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, float> &);
template void EquilibriumDistributionGrid<tags::CPU::Itanium, lbm_applications::LABSWE>::value<double>(double, double,
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, double> &);
