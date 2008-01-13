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

        SPEFrameworkInstruction<1, float, cell::rtm_dma> instruction(oc_element_inverse_float, a.elements(), a.rows() * a.columns());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (MutableMatrix<float>::ElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            if (fabs(*i) <= std::numeric_limits<float>::epsilon())
                *i = 0.0f;
            else
                *i = 1 / *i;
        }

        for (MutableMatrix<float>::ElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            if (fabs(*i) <= std::numeric_limits<float>::epsilon())
                *i = 0.0f;
            else
                *i = 1 / *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }

    DenseVectorContinuousBase<float> &
    ElementInverse<tags::Cell>::value(DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When inverting DenseVectorContinuousBase<float> (Cell):");

        SPEFrameworkInstruction<1, float, cell::rtm_dma> instruction(oc_element_inverse_float, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (Vector<float>::ElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            if (fabs(*i) <= std::numeric_limits<float>::epsilon())
                *i = 0.0f;
            else
                *i = 1 / *i;
        }

        for (Vector<float>::ElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            if (fabs(*i) <= std::numeric_limits<float>::epsilon())
                *i = 0.0f;
            else
                *i = 1 / *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }
}
