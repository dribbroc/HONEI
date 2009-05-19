/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/ppe/spe_instruction.hh>
#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/stringify.hh>

namespace honei
{
    using namespace cell;

    void EquilibriumDistributionGrid<tags::Cell, lbm_applications::LABSWE>::value(float g, float e,
            PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data)
    {
        CONTEXT("When computing LABSWE local equilibrium distribution function (Cell):");

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

        //unsigned long begin((*info.limits)[0]);
        //unsigned long end((*info.limits)[info.limits->size() - 1]);
        /// \todo use accurate begin and end when operation 4 dma supports unaligned data
        unsigned long begin(0);
        unsigned long end(data.h->size());
        unsigned long size(end - begin);

        SPEInstructionQueue q0;
        SPEInstructionQueue q1;
        SPEInstructionQueue q2;
        SPEInstructionQueue q3;
        SPEInstructionQueue q4;

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_0(
                oc_eq_dist_grid_dir_0_float, data.f_eq_0->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e);

        if (instruction_dir_0.use_spe())
        {
            q0.push_back(instruction_dir_0);
        }
        SPEManager::instance()->dispatch(q0);

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_1(
                oc_eq_dist_grid_dir_odd_float, data.f_eq_1->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[1], (*data.distribution_y)[1]);

        if (instruction_dir_1.use_spe())
        {
            q1.push_back(instruction_dir_1);
        }

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_3(
                oc_eq_dist_grid_dir_odd_float, data.f_eq_3->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[3], (*data.distribution_y)[3]);

        if (instruction_dir_3.use_spe())
        {
            q1.push_back(instruction_dir_3);
        }
        SPEManager::instance()->dispatch(q1);

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_5(
                oc_eq_dist_grid_dir_odd_float, data.f_eq_5->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[5], (*data.distribution_y)[5]);

        if (instruction_dir_5.use_spe())
        {
            q2.push_back(instruction_dir_5);
        }

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_7(
                oc_eq_dist_grid_dir_odd_float, data.f_eq_7->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[7], (*data.distribution_y)[7]);

        if (instruction_dir_7.use_spe())
        {
            q2.push_back(instruction_dir_7);
        }
        SPEManager::instance()->dispatch(q2);

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_2(
                oc_eq_dist_grid_dir_even_float, data.f_eq_2->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[2], (*data.distribution_y)[2]);

        if (instruction_dir_2.use_spe())
        {
            q3.push_back(instruction_dir_2);
        }

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_4(
                oc_eq_dist_grid_dir_even_float, data.f_eq_4->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[4], (*data.distribution_y)[4]);

        if (instruction_dir_4.use_spe())
        {
            q3.push_back(instruction_dir_4);
        }
        SPEManager::instance()->dispatch(q3);

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_6(
                oc_eq_dist_grid_dir_even_float, data.f_eq_6->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[6], (*data.distribution_y)[6]);

        if (instruction_dir_6.use_spe())
        {
            q4.push_back(instruction_dir_6);
        }

        SPEFrameworkInstruction<4, float, rtm_dma> instruction_dir_8(
                oc_eq_dist_grid_dir_even_float, data.f_eq_8->elements() + begin, data.h->elements() + begin,
                data.u->elements() + begin, data.v->elements() + begin,
                size, g, e, (*data.distribution_x)[8], (*data.distribution_y)[8]);

        if (instruction_dir_8.use_spe())
        {
            q4.push_back(instruction_dir_8);
        }
            SPEManager::instance()->dispatch(q4);

        // Calculate the last elements on PPU (if needed).
        float e2(e);
        float e23(float(3.) * e2);
        float e26(float(6.) * e2);
        float e42(float(2.) * e2 * e2);
        float e48(float(8.) * e2 * e2);
        float e212(float(12.) * e2);
        float e224(float(24.) * e2);
        for (unsigned long index(begin + instruction_dir_8.transfer_end()) ; index < end ; ++index)
        {
            float u2((*data.u)[index] * (*data.u)[index]);
            float v2((*data.v)[index] * (*data.v)[index]);
            float gh(g * (*data.h)[index]);

            float t1, t2, t3, t4;
            float dxu, dyv;

            // dir 0
            t1 = (float(5.) * gh) / e26;
            t2 = float(2.) / e23 * (u2 + v2);
            (*data.f_eq_0)[index] = (*data.h)[index] * (float(1.) - t1 - t2);

            // dir 1
            dxu = (*data.distribution_x)[1] * (*data.u)[index];
            dyv = (*data.distribution_y)[1] * (*data.v)[index];
            t1 = gh / e26;
            t2 = (dxu + dyv) / e23;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
            t4 = (u2 + v2) / e26;
            (*data.f_eq_1)[index] = (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 3
            dxu = (*data.distribution_x)[3] * (*data.u)[index];
            dyv = (*data.distribution_y)[3] * (*data.v)[index];
            t2 = (dxu + dyv) / e23;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
            (*data.f_eq_3)[index] = (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 5
            dxu = (*data.distribution_x)[5] * (*data.u)[index];
            dyv = (*data.distribution_y)[5] * (*data.v)[index];
            t2 = (dxu + dyv) / e23;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
            (*data.f_eq_5)[index] = (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 7
            dxu = (*data.distribution_x)[7] * (*data.u)[index];
            dyv = (*data.distribution_y)[7] * (*data.v)[index];
            t2 = (dxu + dyv) / e23;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
            (*data.f_eq_7)[index] = (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 2
            dxu = (*data.distribution_x)[2] * (*data.u)[index];
            dyv = (*data.distribution_y)[2] * (*data.v)[index];
            t1 = gh / e224;
            t2 = (dxu + dyv) / e212;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
            t4 = (u2 + v2) / e224;
            (*data.f_eq_2)[index] =  (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 4
            dxu = (*data.distribution_x)[4] * (*data.u)[index];
            dyv = (*data.distribution_y)[4] * (*data.v)[index];
            t2 = (dxu + dyv) / e212;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
            (*data.f_eq_4)[index] =  (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 6
            dxu = (*data.distribution_x)[6] * (*data.u)[index];
            dyv = (*data.distribution_y)[6] * (*data.v)[index];
            t2 = (dxu + dyv) / e212;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
            (*data.f_eq_6)[index] =  (*data.h)[index] * (t1 + t2 + t3 - t4);

            // dir 8
            dxu = (*data.distribution_x)[8] * (*data.u)[index];
            dyv = (*data.distribution_y)[8] * (*data.v)[index];
            t2 = (dxu + dyv) / e212;
            t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
            (*data.f_eq_8)[index] =  (*data.h)[index] * (t1 + t2 + t3 - t4);
        }



        q0.wait();
        q1.wait();
        q2.wait();
        q3.wait();
        q4.wait();

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
}
