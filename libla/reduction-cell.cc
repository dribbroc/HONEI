/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Till Barz <till.barz@uni-dortmund.de>
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
#include <libla/reduction.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

namespace honei
{
    using namespace cell;

    float
    Reduction<rt_sum, tags::Cell>::value(const DenseVector<float> & a)
    {
        CONTEXT("When reducing DenseVector<float> to Scalar by sum (Cell):");

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

        SPEInstruction instruction(oc_dense_float_reduction_sum, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float rest_result(0.0f);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            rest_result += *i;
        }

        if (use_spe)
            instruction.wait();

        return result + rest_result;
    }


    float
    Reduction<rt_min, tags::Cell>::value(const DenseVector<float> & a)
    {
        CONTEXT("When reducing DenseVector<float> to Scalar by minimum (Cell):");

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

        SPEInstruction instruction(oc_dense_float_reduction_min, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }


        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            result = (result > *i)?*i:result;
        }

        if (use_spe)
            instruction.wait();

        return result;
    }

    float
    Reduction<rt_max, tags::Cell>::value(const DenseVector<float> & a)
    {
        CONTEXT("When reducing DenseVector<float> to Scalar by maximum (Cell):");

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

        SPEInstruction instruction(oc_dense_float_reduction_max, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }


        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            result = (result < *i)?*i:result;
        }

        if (use_spe)
            instruction.wait();

        return result;
    }




}
