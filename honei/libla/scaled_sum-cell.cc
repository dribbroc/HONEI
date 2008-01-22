/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/libla/scaled_sum.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/stringify.hh>

namespace honei
{
    using namespace cell;

    DenseVectorContinuousBase<float> &
    ScaledSum<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, const float & c)
    {
        CONTEXT("When adding DenseVectorContinuousBase<float> to scaled DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc, od, oe, of;
        oc.f = c;

        unsigned a_offset((oa.u & 0xF) / sizeof(float)); // Alignment offset of a -> elements calculated on PPU.

        ob.u += (4 * ((4 - a_offset) % 4)); // Adjust SPU start for b respecting the elements calculated on PPU.

        unsigned b_offset((ob.u & 0xF) / sizeof(float)); // Alignment offset of b -> shuffle-factor on SPU.
        of.u = b_offset;

        // Align the address for SPU.
        if (a_offset > 0)
            oa.u += 16 - (4 * a_offset);

        ob.u -= (4 * b_offset);

        od.u = (b.size() - ((4 - a_offset) % 4)) / (1024 * 4); // Subtract PPU-calculated offset from size.
        oe.u = (b.size() - ((4 - a_offset) % 4)) % (1024 * 4);
        oe.u &= ~0xF;

        unsigned rest_index(od.u * 4096 + oe.u + ((4 - a_offset) % 4)); // Rest index for PPU dependent on offset and SPU part.

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

        SPEInstruction instruction(oc_scaled_sum_dense_dense_float, 16 * 1024, oa, ob, oc, od, oe, of);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }
        // Calculate the first a_offset elements on PPU.
        Vector<float>::ConstElementIterator j(b.begin_elements());
        for (Vector<float>::ElementIterator i(a.begin_elements()),
            i_end(a.element_at((4 - a_offset) % 4)) ; i != i_end ; ++i, ++j)
        {
            *i += *j * c;
        }
        // Calculate the last not-by-SPU-done part on PPU.
        Vector<float>::ConstElementIterator k(b.element_at(rest_index));
        Vector<float>::ElementIterator i(a.element_at(rest_index)), i_end(a.end_elements());
        for ( ; i != i_end ; ++i, ++k)
        {
            *i += *k * c;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }

    DenseMatrix<float> &
    ScaledSum<tags::Cell>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b, const float & c)
    {
        CONTEXT("When adding DenseMatrix<float> to scaled DenseMatrix<float> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc, od, oe;
        oc.f = c;
        od.u = (a.rows() * a.columns()) / (1024 * 4);
        oe.u = (a.rows() * a.columns()) % (1024 * 4);
        oe.u &= ~0xF;

        unsigned rest_index(od.u * 4096 + oe.u);

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

        SPEInstruction instruction(oc_scaled_sum_dense_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        Matrix<float>::ConstElementIterator j(b.element_at(rest_index));
        MutableMatrix<float>::ElementIterator i(a.element_at(rest_index)), i_end(a.end_elements());
        for ( ; i != i_end ; ++i, ++j)
        {
            *i += *j * c;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }
}