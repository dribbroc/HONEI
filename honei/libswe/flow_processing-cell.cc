/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/cell/cell.hh>
#include <honei/libswe/flow_processing.hh>
#include <honei/util/memory_backend_cell.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_manager.hh>
#include <honei/util/stringify.hh>

namespace honei
{
    using namespace cell;

    DenseVector<float> &
    FlowProcessing<directions::X, tags::Cell>::value(DenseVector<float> & a)
    {
        CONTEXT("When flow-processing DenseVector<float> for X direction (Cell):");

        SPEFrameworkInstruction<1, float, cell::rtm_dma> instruction(oc_flow_processing_x_float, a.elements(),
                a.size(), 0.0f, 12);

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        if (instruction.use_spe())
            instruction.wait();

        for (Vector<float>::ElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()),
                j(a.element_at(instruction.transfer_end())) ; i.index() + 2 < i_end.index() ; ++i, ++j)
        {
            float h(*i);
            ++i;
            float q1(*i);
            ++i;
            float q2(*i);

            *j = q1;
            ++j;

            *j = q1 * q1 / h + 9.81f * h * h / 2.0f;
            ++j;

            *j = q1 * q2 / h;
        }

        return a;
    }

    DenseVector<float> &
    FlowProcessing<directions::Y, tags::Cell>::value(DenseVector<float> & a)
    {
        CONTEXT("When flow-processing DenseVector<float> for Y direction (Cell):");

        SPEFrameworkInstruction<1, float, cell::rtm_dma> instruction(oc_flow_processing_y_float, a.elements(),
                a.size(), 0.0f, 3);

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        if (instruction.use_spe())
            instruction.wait();

        for (Vector<float>::ElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()),
                j(a.element_at(instruction.transfer_end())) ; i.index() + 2 < i_end.index() ; ++i, ++j)
        {
            float h(*i);
            ++i;
            float q1(*i);
            ++i;
            float q2(*i);

            *j = q2;
            ++j;

            *j = q1 * q2 / h;
            ++j;

            *j = q2 * q2 / h + 9.81f * h * h / 2.0f;
        }

        return a;
    }
}
