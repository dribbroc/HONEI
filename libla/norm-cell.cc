/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
        CONTEXT("When applying max norm to a DenseVector<float> (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od, oe;
        oc.u = a.size() / 4096;
        od.u = a.size() % 4096;
        od.u &= ~0xF;
        oe.f = a[0];
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

        SPEInstruction instruction(oc_dense_float_norm_max, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(a[0]);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (fabs(*i) > fabs(ppu_result)) ? fabs(*i) : fabs(ppu_result);
        }

        if (use_spe)
        {
            instruction.wait();
            result = ppu_result > fabs(result) ? fabs(ppu_result) : fabs(result);
        }
        else
        {
            result = ppu_result;
        }

        return result;
    }


    float
    Norm<vnt_l_two, true, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm to a DenseVector<float> (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od;
        oc.u = a.size() / 4096;
        od.u = a.size() % 4096;
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

        SPEInstruction instruction(oc_dense_float_norm_l_two, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i * *i;
        }

        if (use_spe)
        {
            instruction.wait();
            result = sqrt(ppu_result + result);
        }
        else
        {
            result = sqrt(ppu_result);
        }

        return result;
    }

    float
    Norm<vnt_l_two, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm (quasi-dot-product-version) to a DenseVector<float> (Cell):");
        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od;
        oc.u = a.size() / 4096;
        od.u = a.size() % 4096;
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

        SPEInstruction instruction(oc_dense_float_norm_l_two, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i * *i;
        }

        if (use_spe)
        {
            instruction.wait();
            result = (ppu_result + result);
        }
        else
        {
            result = (ppu_result);
        }

        return result;
    }
}
