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

#include <honei/lbm/extraction_grid.hh>

#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/backends/cell/ppe/spe_instruction.hh>
#include <honei/backends/cell/interface.hh>

#include <honei/la/algorithm.hh>
#include <iostream>

using namespace honei;
using namespace cell;

void ExtractionGrid<tags::Cell, lbm_modes::WET>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon)
{
    CONTEXT("When extracting h, u and v(Cell):");

    //set f to t_temp
    DenseVector<float> * swap;
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

    Operand * operands;
    if (0 != posix_memalign(reinterpret_cast<void **>(&operands), 128, 32 * sizeof(Operand)));

    operands[0].ea = data.h->elements();
    operands[1].ea = data.u->elements();
    operands[2].ea = data.v->elements();
    operands[3].ea = data.f_0->elements();
    operands[4].ea = data.f_1->elements();
    operands[5].ea = data.f_2->elements();
    operands[6].ea = data.f_3->elements();
    operands[7].ea = data.f_4->elements();
    operands[8].ea = data.f_5->elements();
    operands[9].ea = data.f_6->elements();
    operands[10].ea = data.f_7->elements();
    operands[11].ea = data.f_8->elements();
    operands[12].fa[0] = (*data.distribution_x)[0];
    operands[13].fa[0] = (*data.distribution_x)[1];
    operands[14].fa[0] = (*data.distribution_x)[2];
    operands[15].fa[0] = (*data.distribution_x)[3];
    operands[16].fa[0] = (*data.distribution_x)[4];
    operands[17].fa[0] = (*data.distribution_x)[5];
    operands[18].fa[0] = (*data.distribution_x)[6];
    operands[19].fa[0] = (*data.distribution_x)[7];
    operands[20].fa[0] = (*data.distribution_x)[8];
    operands[21].fa[1] = (*data.distribution_y)[0];
    operands[22].fa[1] = (*data.distribution_y)[1];
    operands[23].fa[1] = (*data.distribution_y)[2];
    operands[24].fa[1] = (*data.distribution_y)[3];
    operands[25].fa[1] = (*data.distribution_y)[4];
    operands[26].fa[1] = (*data.distribution_y)[5];
    operands[27].fa[1] = (*data.distribution_y)[6];
    operands[28].fa[1] = (*data.distribution_y)[7];
    operands[29].fa[1] = (*data.distribution_y)[8];

    unsigned long begin(0);
    unsigned long end(data.h->size());
    unsigned long size(end - begin);

    SPEFrameworkInstruction<1, float, rtm_dma> instruction(
            oc_extraction_grid_float, data.f_0->elements(), size);
    instruction.instruction().a.ea = operands;

    SPEManager::instance()->dispatch(instruction);

    for (unsigned long i(begin + instruction.transfer_end()) ; i < end ; ++i)
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
    instruction.wait();
}
