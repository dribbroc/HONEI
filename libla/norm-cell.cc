/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007,2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Till Barz <till.barz@uni-dortmund.de>
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
#include <libla/norm.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

#include <cmath>

namespace honei
{
    using namespace cell;

    float
    Norm<vnt_max, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying max norm to DenseVectorContinuousBase<float> (Cell):");

        float result(0.0f);

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_norm_max_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result > fabs(*i)) ? ppu_result : fabs(*i);
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result > fabs(*i)) ? ppu_result : fabs(*i);
        }

        if (instruction.use_spe())
            instruction.wait();

        return (ppu_result > result) ? ppu_result : result;
    }

    float
    Norm<vnt_l_one, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L1-norm to DenseVectorContinuousBase<float> (Cell):");

        float result(0.0f);

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_norm_l_one_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            ppu_result += fabs(*i);
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result += fabs(*i);
        }

        if (instruction.use_spe())
            instruction.wait();

        return result + ppu_result;
    }

    float
    Norm<vnt_l_two, true, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm (with square-root) to DenseVectorContinuousBase<float> (Cell):");

        return sqrt(Norm<vnt_l_two, false, tags::Cell>::value(a));
    }

    float
    Norm<vnt_l_two, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm to DenseVectorContinuousBase<float> (Cell):");

        float result(0.0f);

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_norm_l_two_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            float temp(fabs(*i));
            temp *= temp;

            ppu_result += temp;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            float temp(fabs(*i));
            temp *= temp;

            ppu_result += temp;
        }

        if (instruction.use_spe())
            instruction.wait();

        return result + ppu_result;
    }
}

