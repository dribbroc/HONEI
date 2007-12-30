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
#include <libla/element_inverse.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float> &
    ElementInverse<tags::Cell>::value(DenseMatrix<float> & a)
    {
        CONTEXT("When inverting DenseMatrix<float> (Cell):");

        Operand oa = { a.elements() };
        Operand ob, oc;
        ob.u = (a.rows() * a.columns()) / 4096;
        oc.u = (a.rows() * a.columns()) % 4096;
        oc.u &= ~0xF;

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

        SPEInstruction instruction(oc_element_inverse_float, 16 * 1024, oa, ob, oc);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (MutableMatrix<float>::ElementIterator j(a.element_at(rest_index)), j_end(a.end_elements()) ; j != j_end ; ++j)
        {
            if (*j == 0)
                continue;

            *j = 1 / *j;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }

    DenseVectorContinuousBase<float> &
    ElementInverse<tags::Cell>::value(DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When inverting DenseVectorContinuousBase<float> (Cell):");

        Operand oa = { a.elements() };
        Operand ob, oc;

        unsigned offset((oa.u & 0xF) / sizeof(float));
        if (a.size() < 5)
            offset = 0;

        if (offset > 0)
            oa.u += 16 -(4 * offset); // Align the address for SPU.


        ob.u = (a.size() - ((4 - offset) % 4)) / (1024 * 4); // Subtract PPU-calculated offset from size.
        oc.u = (a.size() - ((4 - offset) % 4)) % (1024 * 4);
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

        SPEInstruction instruction(oc_element_inverse_float, 16 * 1024, oa, ob, oc);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (Vector<float>::ElementIterator i(a.begin_elements()),
                i_end(a.element_at((4 - offset) % 4)) ; i != i_end ; ++i)
        {
            if (*i == 0.0f)
                continue;

            *i = 1 / *i;
        }

        for (Vector<float>::ElementIterator j(a.element_at(rest_index)), j_end(a.end_elements()) ; j != j_end ; ++j)
        {
            if (*j == 0.0f)
                continue;

            *j = 1 / *j;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }


}
