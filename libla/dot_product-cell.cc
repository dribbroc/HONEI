/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007. 2008 Sven Mallach <sven.mallach@honei.org>
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

#include <cell/cell.hh>
#include <libla/dot_product.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace honei
{
    using namespace cell;

    float
    DotProduct<tags::Cell>::value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When adding DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        float result(0.0f);
        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc = { b.elements() };
        Operand od, oe, of, og, oh, oi, oj;

        unsigned a_offset((ob.u & 0xF) / sizeof(float)); // Alignment offset of a -> elements calculated on PPU.
        if (a.size() < 5)
            a_offset = 0;


        unsigned skip = ((4 - a_offset) % 4);
        oc.u += (4 * skip); // Adjust SPU start for b respecting the elements calculated on PPU.

        unsigned b_offset((oc.u & 0xF) / sizeof(float)); // Alignment offset of b -> shuffle-factor on SPU.
        of.u = b_offset;

        // Align the address for SPU.
        ob.u += (4 * skip);
        oc.u += (4 * (4 - b_offset));

        float * dma_start = reinterpret_cast<float *>(oc.u);
        og.f = *(dma_start - 4);
        oh.f = *(dma_start - 3);
        oi.f = *(dma_start - 2);
        oj.f = *(dma_start - 1);

        //Subtract PPU-calculated parts from size.
        od.u = (b.size() - skip) / (1024 * 4);
        oe.u = (b.size() - skip) % (1024 * 4);
        oe.u &= ~0xF;

        // Rest index dependent on offset and SPU part.
        unsigned rest_index(od.u * 4096 + oe.u + skip);

        oe.u *= 4;

        bool use_spe(true);

        if (0 == oe.u)
        {
            if (od.u > 0)
            {
                oe.u = 16 * 1024;
            }
            else
            {
                use_spe = false;
            }
        }
        else
        {
            ++od.u;
        }
        SPEInstruction instruction(oc_dot_product_dense_dense_float, 16 * 1024, oa, ob, oc, od, oe, of, og, oh, oi, oj);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float rest_result(0.0f);

        // Calculate the first 4 - a_offset elements on PPU.
        Vector<float>::ConstElementIterator j(b.begin_elements());
        for (Vector<float>::ConstElementIterator i(a.begin_elements()),
            i_end(a.element_at(skip)) ; i != i_end ; ++i, ++j)
        {
            rest_result += *i * *j;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()),
                j(b.element_at(rest_index)) ; i != i_end ; ++i, ++j)
        {
            rest_result += *i * *j;
        }

        if (use_spe)
            instruction.wait();

        return result + rest_result;
    }
}
