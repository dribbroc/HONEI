/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
    float
    DotProduct<tags::Cell>::value(const DenseVector<float> & a, const DenseVector<float> & b)
    {
        CONTEXT("When adding DenseVector<float> to DenseVector<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        float result(0.0f);
        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc = { b.elements() };
        Operand od, oe;
        // hardcode transfer buffer size for now.
        od.u = a.size() / (1024 * 4);
        oe.u = a.size() % (1024 * 4);
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

        SPEInstruction instruction(oc_dense_dense_float_dot_product, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float rest_result(0.0f);

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
