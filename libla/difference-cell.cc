/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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
#include <libla/difference.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

namespace honei
{
    using namespace cell;

    DenseVectorContinuousBase<float> &
    Difference<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When subtracting DenseVectorContinuousBase<float> from DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc, od, oe;

        unsigned a_offset((oa.u & 0xF) / sizeof(float)); // Alignment offset of a -> elements calculated on PPU.
        if (a.size() < 5)
            a_offset = 0;

        ob.u += (4 * ((4 - a_offset) % 4)); // Adjust SPU start for b respecting the elements calculated on PPU.

        unsigned b_offset((ob.u & 0xF) / sizeof(float)); // Alignment offset of b -> shuffle-factor on SPU.
        oe.u = b_offset;

        // Align the address for SPU.
        if (a_offset > 0)
            oa.u += 16 - (4 * a_offset);

        ob.u -= (4 * b_offset);

        oc.u = (b.size() - ((4 - a_offset) % 4)) / (1024 * 4); // Subtract PPU-calculated offset from size.
        od.u = (b.size() - ((4 - a_offset) % 4)) % (1024 * 4);
        od.u &= ~0xF;

        unsigned rest_index(oc.u * 4096 + od.u + ((4 - a_offset) % 4)); // Rest index dependent on offset and SPU part.

        od.u *= 4;

        bool use_spe(true);

        if (0 == od.u)
        {
            if (oc.u > 0)
            {
                od.u = 16 * 1024;
            }
            else
            {
                use_spe = false;
            }
        }
        else
        {
            ++oc.u;
        }

        SPEInstruction instruction(oc_difference_dense_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        // Calculate the first a_offset elements on PPU.
        Vector<float>::ConstElementIterator j(b.begin_elements());
        for (Vector<float>::ElementIterator i(a.begin_elements()),
                i_end(a.element_at((4 - a_offset) % 4)) ; i != i_end ; ++i, ++j)
        {
            *i -= *j;
        }

        Vector<float>::ConstElementIterator k(b.element_at(rest_index));
        for (Vector<float>::ElementIterator i(a.element_at(rest_index)),
                i_end(a.end_elements()) ; i != i_end ; ++i, ++k)
        {
            *i -= *k;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }

    DenseVector<float> &
    Difference<tags::Cell>::value(const SparseVector<float> & a, DenseVector<float> & b)
    {
        CONTEXT("When subtracting SparseVector<float> from DenseVector<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        Scale<tags::Cell>::value(-1.0f, b);
        Sum<tags::Cell>::value(b, a);

        return b;
    }

}
