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
#include <libla/sum.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

namespace honei
{
    DenseVector<float> &
    Sum<tags::Cell>::value(DenseVector<float> & a, const DenseVector<float> & b)
    {
        CONTEXT("When adding DenseVector<float> to DenseVector<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc, od;
        // hardcode transfer buffer size for now.
        oc.u = a.size() / (1024 * 4);
        od.u = a.size() % (1024 * 4);
        od.u &= ~0xF;

        unsigned rest_index(oc.u * 4096 + od.u);

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

        SPEInstruction instruction(oc_dense_dense_float_sum, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        Vector<float>::ConstElementIterator j(b.element_at(rest_index));
        Vector<float>::ElementIterator i(a.element_at(rest_index)), i_end(a.end_elements());
        for ( ; i != i_end ; ++i, ++j)
        {
            *i += *j;
        }

        if (use_spe)
            instruction.wait();

        return a;
    }

    DenseVector<float> &
    Sum<tags::Cell>::value(DenseVector<float> & a, const SparseVector<float> & b)
    {
        CONTEXT("When adding DenseVector<float> to SparseVector<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc = { b.indices() };
        Operand od;
        od.u = b.used_elements();
        SPEInstruction instruction(oc_dense_sparse_float_sum, a.size(), oa, ob, oc, od);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();

        return a;
    }


}
