/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/cell/ppe/spe_transfer_list.hh>

#include <honei/lbm/collide_stream_grid.hh>

#include <iostream>


using namespace honei;

void CollideStreamGrid<tags::Cell, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                float tau)
{
    using namespace cell;
    CONTEXT("When performing collision and streaming (Cell):");

    unsigned long * limits(info.limits->elements());
    unsigned long * dir_index_1(info.dir_index_1->elements());
    unsigned long * dir_index_2(info.dir_index_2->elements());
    unsigned long * dir_index_3(info.dir_index_3->elements());
    unsigned long * dir_index_4(info.dir_index_4->elements());
    unsigned long * dir_index_5(info.dir_index_5->elements());
    unsigned long * dir_index_6(info.dir_index_6->elements());
    unsigned long * dir_index_7(info.dir_index_7->elements());
    unsigned long * dir_index_8(info.dir_index_8->elements());
    unsigned long * dir_1(info.dir_1->elements());
    unsigned long * dir_2(info.dir_2->elements());
    unsigned long * dir_3(info.dir_3->elements());
    unsigned long * dir_4(info.dir_4->elements());
    unsigned long * dir_5(info.dir_5->elements());
    unsigned long * dir_6(info.dir_6->elements());
    unsigned long * dir_7(info.dir_7->elements());
    unsigned long * dir_8(info.dir_8->elements());
    float * f_temp_0(data.f_temp_0->elements());
    float * f_temp_1(data.f_temp_1->elements());
    float * f_temp_2(data.f_temp_2->elements());
    float * f_temp_3(data.f_temp_3->elements());
    float * f_temp_4(data.f_temp_4->elements());
    float * f_temp_5(data.f_temp_5->elements());
    float * f_temp_6(data.f_temp_6->elements());
    float * f_temp_7(data.f_temp_7->elements());
    float * f_temp_8(data.f_temp_8->elements());
    float * f_0(data.f_0->elements());
    float * f_1(data.f_1->elements());
    float * f_2(data.f_2->elements());
    float * f_3(data.f_3->elements());
    float * f_4(data.f_4->elements());
    float * f_5(data.f_5->elements());
    float * f_6(data.f_6->elements());
    float * f_7(data.f_7->elements());
    float * f_8(data.f_8->elements());
    float * f_eq_0(data.f_eq_0->elements());
    float * f_eq_1(data.f_eq_1->elements());
    float * f_eq_2(data.f_eq_2->elements());
    float * f_eq_3(data.f_eq_3->elements());
    float * f_eq_4(data.f_eq_4->elements());
    float * f_eq_5(data.f_eq_5->elements());
    float * f_eq_6(data.f_eq_6->elements());
    float * f_eq_7(data.f_eq_7->elements());
    float * f_eq_8(data.f_eq_8->elements());

    SPEInstructionQueue q0;
    SPEInstructionQueue q1;
    SPEInstructionQueue q2;
    SPEInstructionQueue q3;
    SPEInstructionQueue q4;
    SPEInstructionQueue q5;
    unsigned a_offset, skip, width;
    width = dir_1[0] + dir_index_1[1];

    a_offset = (((unsigned long)f_temp_1 + dir_1[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_1[0];
    skip+= dir_index_1[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_1(oc_collide_stream_grid_times_float,
            f_temp_1 + dir_1[0], f_1 + dir_index_1[0], f_eq_1 + dir_index_1[0], dir_index_1[info.dir_index_1->size() - 1] - dir_index_1[0], tau, skip, width);
    if (instruction_1.use_spe())
    {
        q1.push_back(instruction_1);
    }

    a_offset = (((unsigned long)f_temp_2 + dir_2[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_2[0];
    skip+= dir_index_2[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_2(oc_collide_stream_grid_times_float,
            f_temp_2 + dir_2[0], f_2 + dir_index_2[0], f_eq_2 + dir_index_2[0], dir_index_2[info.dir_index_2->size() - 1] - dir_index_2[0], tau, skip, width);
    if (instruction_2.use_spe())
    {
        q1.push_back(instruction_2);
    }
    SPEManager::instance()->dispatch(q1);


    a_offset = (((unsigned long)f_temp_4 + dir_4[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_4[0];
    skip+= dir_index_4[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_4(oc_collide_stream_grid_times_float,
            f_temp_4 + dir_4[0], f_4 + dir_index_4[0], f_eq_4 + dir_index_4[0], dir_index_4[info.dir_index_4->size() - 1] - dir_index_4[0], tau, skip, width);
    if (instruction_4.use_spe())
    {
        q2.push_back(instruction_4);
    }
    SPEManager::instance()->dispatch(q2);

    a_offset = (((unsigned long)f_temp_5 + dir_5[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_5[0];
    skip+= dir_index_5[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_5(oc_collide_stream_grid_times_float,
            f_temp_5 + dir_5[0], f_5 + dir_index_5[0], f_eq_5 + dir_index_5[0], dir_index_5[info.dir_index_5->size() - 1] - dir_index_5[0], tau, skip, width);
    if (instruction_5.use_spe())
    {
        q3.push_back(instruction_5);
    }
    SPEManager::instance()->dispatch(q3);


    a_offset = (((unsigned long)f_temp_6 + dir_6[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_6[0];
    skip+= dir_index_6[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_6(oc_collide_stream_grid_times_float,
            f_temp_6 + dir_6[0], f_6 + dir_index_6[0], f_eq_6 + dir_index_6[0], dir_index_6[info.dir_index_6->size() - 1] - dir_index_6[0], tau, skip, width);
    if (instruction_6.use_spe())
    {
        q4.push_back(instruction_6);
    }
    SPEManager::instance()->dispatch(q4);


    a_offset = (((unsigned long)f_temp_8 + dir_8[0] * sizeof(float))& 0xF) / sizeof(float);
    skip = (16 / sizeof(float) - a_offset) % (16 / sizeof(float));
    skip+= dir_8[0];
    skip+= dir_index_8[0];
    SPEFrameworkInstruction<3, float, rtm_dma> instruction_8(oc_collide_stream_grid_times_float,
            f_temp_8 + dir_8[0], f_8 + dir_index_8[0], f_eq_8 + dir_index_8[0], dir_index_8[info.dir_index_8->size() - 1] - dir_index_8[0], tau, skip, width);
    if (instruction_8.use_spe())
    {
        q5.push_back(instruction_8);
    }
    SPEManager::instance()->dispatch(q5);

    SPEFrameworkInstruction<3, float, rtm_dma> instruction_3(oc_collide_stream_grid_float,
            f_temp_3 + dir_3[0], f_3 + dir_index_3[0], f_eq_3 + dir_index_3[0], dir_index_3[info.dir_index_3->size() - 1] - dir_index_3[0], tau);
    if (instruction_3.use_spe())
    {
        q0.push_back(instruction_3);
    }

    SPEFrameworkInstruction<3, float, rtm_dma> instruction_7(oc_collide_stream_grid_float,
            f_temp_7 + dir_7[0], f_7 + dir_index_7[0], f_eq_7 + dir_index_7[0], dir_index_7[info.dir_index_7->size() - 1] - dir_index_7[0], tau);
    if (instruction_7.use_spe())
    {
        q0.push_back(instruction_7);
    }

    SPEFrameworkInstruction<3, float, rtm_dma> instruction_0(oc_collide_stream_grid_float,
            f_temp_0, f_0, f_eq_0, limits[info.limits->size() - 1], tau);
    if (instruction_0.use_spe())
    {
        q0.push_back(instruction_0);
    }
    SPEManager::instance()->dispatch(q0);

    for (unsigned long i(instruction_0.transfer_end()) ; i < limits[info.limits->size() - 1] ; ++i)
    {
        f_temp_0[i] = f_0[i] - (f_0[i] - f_eq_0[i])/tau;
    }
    for (unsigned long i(dir_index_1[0]), offset(dir_1[0]) ; offset < dir_1[0] + instruction_1.transfer_begin() ; ++i, ++offset)
    {
        f_temp_1[offset] = f_1[i] - (f_1[i] - f_eq_1[i])/tau;
    }
    for (unsigned long i(dir_index_1[0] + instruction_1.transfer_end()), offset(dir_1[0] + instruction_1.transfer_end()) ; i < dir_index_1[info.dir_index_1->size() - 1] ; ++i, ++offset)
    {
        f_temp_1[offset] = f_1[i] - (f_1[i] - f_eq_1[i])/tau;
    }
    for (unsigned long i(dir_index_2[0]), offset(dir_2[0]) ; offset < dir_2[0] + instruction_2.transfer_begin() ; ++i, ++offset)
    {
        f_temp_2[offset] = f_2[i] - (f_2[i] - f_eq_2[i])/tau;
    }
    for (unsigned long i(dir_index_2[0] + instruction_2.transfer_end()), offset(dir_2[0] + instruction_2.transfer_end()) ; i < dir_index_2[info.dir_index_2->size() - 1] ; ++i, ++offset)
    {
        f_temp_2[offset] = f_2[i] - (f_2[i] - f_eq_2[i])/tau;
    }
    for (unsigned long i(dir_index_3[0]), offset(dir_3[0]) ; offset < dir_3[0] + instruction_3.transfer_begin() ; ++i, ++offset)
    {
        f_temp_3[offset] = f_3[i] - (f_3[i] - f_eq_3[i])/tau;
    }
    for (unsigned long i(dir_index_3[0] + instruction_3.transfer_end()), offset(dir_3[0] + instruction_3.transfer_end()) ; i < dir_index_3[info.dir_index_3->size() - 1] ; ++i, ++offset)
    {
        f_temp_3[offset] = f_3[i] - (f_3[i] - f_eq_3[i])/tau;
    }
    for (unsigned long i(dir_index_4[0]), offset(dir_4[0]) ; offset < dir_4[0] + instruction_4.transfer_begin() ; ++i, ++offset)
    {
        f_temp_4[offset] = f_4[i] - (f_4[i] - f_eq_4[i])/tau;
    }
    for (unsigned long i(dir_index_4[0] + instruction_4.transfer_end()), offset(dir_4[0] + instruction_4.transfer_end()) ; i < dir_index_4[info.dir_index_4->size() - 1] ; ++i, ++offset)
    {
        f_temp_4[offset] = f_4[i] - (f_4[i] - f_eq_4[i])/tau;
    }
    for (unsigned long i(dir_index_5[0]), offset(dir_5[0]) ; offset < dir_5[0] + instruction_5.transfer_begin() ; ++i, ++offset)
    {
        f_temp_5[offset] = f_5[i] - (f_5[i] - f_eq_5[i])/tau;
    }
    for (unsigned long i(dir_index_5[0] + instruction_5.transfer_end()), offset(dir_5[0] + instruction_5.transfer_end()) ; i < dir_index_5[info.dir_index_5->size() - 1] ; ++i, ++offset)
    {
        f_temp_5[offset] = f_5[i] - (f_5[i] - f_eq_5[i])/tau;
    }
    for (unsigned long i(dir_index_6[0]), offset(dir_6[0]) ; offset < dir_6[0] + instruction_6.transfer_begin() ; ++i, ++offset)
    {
        f_temp_6[offset] = f_6[i] - (f_6[i] - f_eq_6[i])/tau;
    }
    for (unsigned long i(dir_index_6[0] + instruction_6.transfer_end()), offset(dir_6[0] + instruction_6.transfer_end()) ; i < dir_index_6[info.dir_index_6->size() - 1] ; ++i, ++offset)
    {
        f_temp_6[offset] = f_6[i] - (f_6[i] - f_eq_6[i])/tau;
    }
    for (unsigned long i(dir_index_7[0]), offset(dir_7[0]) ; offset < dir_7[0] + instruction_7.transfer_begin() ; ++i, ++offset)
    {
        f_temp_7[offset] = f_7[i] - (f_7[i] - f_eq_7[i])/tau;
    }
    for (unsigned long i(dir_index_7[0] + instruction_7.transfer_end()), offset(dir_7[0] + instruction_7.transfer_end()) ; i < dir_index_7[info.dir_index_7->size() - 1] ; ++i, ++offset)
    {
        f_temp_7[offset] = f_7[i] - (f_7[i] - f_eq_7[i])/tau;
    }
    for (unsigned long i(dir_index_8[0]), offset(dir_8[0]) ; offset < dir_8[0] + instruction_8.transfer_begin() ; ++i, ++offset)
    {
        f_temp_8[offset] = f_8[i] - (f_8[i] - f_eq_8[i])/tau;
    }
    for (unsigned long i(dir_index_8[0] + instruction_8.transfer_end()), offset(dir_8[0] + instruction_8.transfer_end()) ; i < dir_index_8[info.dir_index_8->size() - 1] ; ++i, ++offset)
    {
        f_temp_8[offset] = f_8[i] - (f_8[i] - f_eq_8[i])/tau;
    }


    q0.wait();
    q1.wait();
    q2.wait();
    q3.wait();
    q4.wait();
    q5.wait();

}
