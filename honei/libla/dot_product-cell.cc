/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/libla/dot_product.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/partitioner.hh>

namespace honei
{
    using namespace cell;

    float
    DotProduct<tags::Cell>::value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When adding DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        float ppu_result(0.0f), spu_result(0.0f);

        SPEFrameworkInstruction<2, float, rtm_mail> instruction(oc_dot_product_dense_dense_float, &spu_result,
                        a.elements(), b.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < instruction.transfer_begin() ; ++index)
        {
            ppu_result += a[index] * b[index];
        }

        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(instruction.transfer_end()) ; index < a.size() ; index++)
        {
            ppu_result += a[index] * b[index];
        }

        if (instruction.use_spe())
        {
            instruction.wait();
        }

        return ppu_result + spu_result;
    }
}
