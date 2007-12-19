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

#include <cell/cell.hh>
#include <libla/scale.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float> &
    Scale<tags::Cell>::value(const float a, DenseMatrix<float> & b)
    {
        CONTEXT("When scaling DenseMatrix<float> (Cell):");

        Operand oa = { b.elements() };
        Operand ob;
        ob.u = (b.rows() * b.columns()) / (1024 * 4);
        Operand oc;
        oc.u = (b.rows() * b.columns()) % (1024 * 4);
        oc.u &= ~0xF;
        Operand od;
        od.f = a;

        unsigned rest_index(ob.u * 4096 + oc.u);

        oc.u *= 4;

        bool use_spe(true);

        if (0 == oc.u)
        {
            if (ob.u > 0)
            {
                oc.u = 16 * 1024;
            }
            else
            {
                use_spe = false;
            }
        }
        else
        {
            ++ob.u;
        }

        SPEInstruction instruction(oc_dense_float_scale, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        MutableMatrix<float>::ElementIterator i(b.element_at(rest_index)), i_end(b.end_elements());
        for ( ; i != i_end ; ++i)
        {
            *i *= od.f;
        }

        if (use_spe)
            instruction.wait();

        return b;
    }

    DenseVectorContinuousBase<float> &
    Scale<tags::Cell>::value(const float a, DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When scaling DenseVector<float> (Cell):");

        Operand oa = { b.elements() };
        Operand ob, oc, od;
        od.f = a;

        unsigned offset((oa.u & 0xF) / sizeof(float));

        if (offset > 0)
            oa.u += 16 -(4 * offset); // Align the address for SPU.

        ob.u = (b.size() - ((4 - offset) % 4)) / (1024 * 4); // Subtract PPU-calculated offset from size.
        oc.u = (b.size() - ((4 - offset) % 4)) % (1024 * 4);
        oc.u &= ~0xF;

        unsigned rest_index(ob.u * 4096 + oc.u + ((4 - offset) % 4)); // Rest index for PPU dependent on offset and SPU part.
        oc.u *= 4;

        bool use_spe(true);

        if (0 == oc.u)
        {
            if (ob.u > 0)
            {
                oc.u = 16 * 1024;
            }
            else
            {
                use_spe = false;
            }
        }
        else
        {
            ++ob.u;
        }

        SPEInstruction instruction(oc_dense_float_scale, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (Vector<float>::ElementIterator i(b.begin_elements()),
                i_end(b.element_at((4 - offset) % 4)) ; i != i_end ; ++i)
        {
            *i *= od.f;
        }

        for (Vector<float>::ElementIterator i(b.element_at(rest_index)),
                i_end(b.end_elements()) ; i != i_end ; ++i)
        {
            *i *= od.f;
        }

        if (use_spe)
            instruction.wait();

        return b;
    }

    SparseVector<float> &
    Scale<tags::Cell>::value(const float a, SparseVector<float> & b)
    {
        CONTEXT("When scaling SparseVector<float> (Cell):");

        Operand oa = { b.elements() };
        Operand ob;
        ob.u = b.size() / (1024 * 4);
        Operand oc;
        oc.u = b.size() % (1024 * 4);
        oc.u &= ~0xF;
        Operand od;
        od.f = a;

        unsigned rest_index(ob.u * 4096 + oc.u);

        oc.u *= 4;

        bool use_spe(true);

        if (0 == oc.u)
        {
            if (ob.u > 0)
            {
                oc.u = 16 * 1024;
            }
            else
            {
                use_spe = false;
            }
        }
        else
        {
            ++ob.u;
        }

        SPEInstruction instruction(oc_dense_float_scale, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (unsigned i(rest_index); i < b.used_elements() ; i++)
        {
            *(b.elements() + i) *= a;
        }

        if (use_spe)
            instruction.wait();

        return b;
    }
}

