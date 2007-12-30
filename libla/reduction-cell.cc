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
#include <libla/reduction.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>

namespace honei
{
    using namespace cell;

    float
    Reduction<rt_sum, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by sum (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od;

        unsigned offset((ob.u & 0xF) / sizeof(float));

        if (offset > 0)
            ob.u += 16 -(4 * offset); // Align the address for SPU.

        if (a.size() < 5)
            offset = 0;

        oc.u = (a.size() - ((4 - offset) % 4)) / (1024 * 4); // Subtract PPU-calculated offset from size.
        od.u = (a.size() - ((4 - offset) % 4)) % (1024 * 4);
        od.u &= ~0xF;

        unsigned rest_index(oc.u * 4096 + od.u + ((4 - offset) % 4)); // Rest index for PPU dependent on offset and SPU part.
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
        SPEInstruction instruction(oc_reduction_sum_dense_float, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at((4 - offset) % 4)) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        if (use_spe)
            instruction.wait();

        return result += ppu_result;
    }

    DenseVector<float>
    Reduction<rt_sum, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by sum (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (unsigned i(0) ; i < a.rows() ; i++)
        {
            Operand oa = { result.elements() + i };
            Operand ob = { a.elements() + (i * a.columns()) };
            Operand oc, od;
            oc.u = a.columns() / 4096;
            od.u = a.columns() % 4096;
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

            SPEInstruction instruction(oc_reduction_sum_dense_float, 16 * 1024, oa, ob, oc, od);

            if (use_spe)
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(0.0f);
            for (Vector<float>::ConstElementIterator j(a[i].element_at(rest_index)), j_end(a[i].end_elements()) ; j != j_end ; ++j)
            {
                ppu_result += *j;
            }

            if (use_spe)
                instruction.wait();

            result[i] += ppu_result;

        }
        return result;
    }

    float
    Reduction<rt_min, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by minimum (Cell):");

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

        SPEInstruction instruction(oc_reduction_min_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
            unsigned offset((ob.u & 0xF) / sizeof(float));
            rest_index -=offset;
        }

        float ppu_result(a[0]);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (*i < ppu_result) ? *i : ppu_result;
        }

        if (use_spe)
        {
            instruction.wait();
            result = ppu_result < result ? ppu_result : result;
        }
        else
        {
            result = ppu_result;
        }

        return result;
    }

    DenseVector<float>
    Reduction<rt_min, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by min (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (unsigned i(0) ; i < a.rows() ; i++)
        {
            Operand oa = { result.elements() + i };
            Operand ob = { a.elements() + (i * a.columns()) };
            Operand oc, od, oe;
            oc.u = a.columns() / 4096;
            od.u = a.columns() % 4096;
            od.u &= ~0xF;
            oe.f = *(a.elements() + (i * a.columns()));

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

            SPEInstruction instruction(oc_reduction_min_dense_float, 16 * 1024, oa, ob, oc, od, oe);

            if (use_spe)
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(*(a.elements() + (i * a.columns())));
            for (Vector<float>::ConstElementIterator j(a[i].element_at(rest_index)), j_end(a[i].end_elements()) ; j != j_end ; ++j)
            {
                ppu_result = (*j < ppu_result) ? *j : ppu_result;
            }

            if (use_spe)
            {
                instruction.wait();
                result[i] = ppu_result < result[i] ? ppu_result : result[i];
            }
            else
            {
               result[i] = ppu_result;
            }

       }
        return result;
    }

    float
    Reduction<rt_max, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by maximum (Cell):");

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

        SPEInstruction instruction(oc_reduction_max_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
            unsigned offset((ob.u & 0xF) / sizeof(float));
            rest_index -=offset;
        }

        float ppu_result(a[0]);
        for (Vector<float>::ConstElementIterator i(a.element_at(rest_index)), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (*i > ppu_result) ? *i : ppu_result;
        }

        if (use_spe)
        {
            instruction.wait();
            result = ppu_result > result ? ppu_result : result;
        }
        else
        {
            result = ppu_result;
        }

        return result;
    }

    DenseVector<float>
    Reduction<rt_max, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by max (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (unsigned i(0) ; i < a.rows() ; i++)
        {
            Operand oa = { result.elements() + i };
            Operand ob = { a.elements() + (i * a.columns()) };
            Operand oc, od, oe;
            oc.u = a.columns() / 4096;
            od.u = a.columns() % 4096;
            od.u &= ~0xF;
            oe.f = *(a.elements() + (i * a.columns()));
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

            SPEInstruction instruction(oc_reduction_max_dense_float, 16 * 1024, oa, ob, oc, od, oe);

            if (use_spe)
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(*(a.elements() + (i * a.columns())));
            for (Vector<float>::ConstElementIterator j(a[i].element_at(rest_index)), j_end(a[i].end_elements()) ; j != j_end ; ++j)
            {
                ppu_result = (*j > ppu_result) ? *j : ppu_result;
            }

            if (use_spe)
            {
                instruction.wait();
                result[i] = ppu_result > result[i] ? ppu_result : result[i];
            }
            else
            {
               result[i] = ppu_result;
            }

        }
        return result;
    }
   float
    Reduction<rt_sum, tags::Cell>::value(const SparseVector<float> & a)
    {
        CONTEXT("When reducing SparseVector<float> to Scalar by sum (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od;
        oc.u = a.used_elements() / 4096;
        od.u = a.used_elements() % 4096;
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

        SPEInstruction instruction(oc_reduction_sum_dense_float, 16 * 1024, oa, ob, oc, od);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);
        for (unsigned i(rest_index); i < a.used_elements() ; i++)
        {
            ppu_result += *(a.elements() + i);
        }

        if (use_spe)
            instruction.wait();

        return (result + ppu_result);
    }

    DenseVector<float>
    Reduction<rt_sum, tags::Cell>::value(const SparseMatrix<float> & a)
    {
        CONTEXT("When reducing SparseMatrix<float> to Vector by sum (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (SparseMatrix<float>::ConstRowIterator i(a.begin_non_zero_rows()),
                i_end(a.end_non_zero_rows()) ; i != i_end ; ++i)
        {
            result[i.index()] = Reduction<rt_sum, tags::Cell>::value(*i);
        }

        return result;
    }

    float
    Reduction<rt_min, tags::Cell>::value(const SparseVector<float> & a)
    {
        CONTEXT("When reducing SparseVector<float> to Scalar by minimum (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od, oe;
        oc.u = a.used_elements() / 4096;
        od.u = a.used_elements() % 4096;
        od.u &= ~0xF;
        oe.f = a.used_elements() < a.size() ? 0.0f : a[*a.indices()];

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

        SPEInstruction instruction(oc_reduction_min_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(oe.f);

        for (unsigned i(rest_index); i < a.used_elements() ; i++)
        {
            ppu_result = *(a.elements() + i) < ppu_result ? *(a.elements() + i) : ppu_result;
        }

        if (use_spe)
        {
            instruction.wait();
            result = ppu_result < result ? ppu_result : result;
        }
        else
        {
            result = ppu_result;
        }

        return result;
    }

    DenseVector<float>
    Reduction<rt_min, tags::Cell>::value(const SparseMatrix<float> & a)
    {
        CONTEXT("When reducing SparseMatrix<float> to Vector by min (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (SparseMatrix<float>::ConstRowIterator i(a.begin_non_zero_rows()), 
                i_end(a.end_non_zero_rows()) ; i != i_end ; ++i)
        {
            result[i.index()] = Reduction<rt_min, tags::Cell>::value(*i);
        }

        return result;
    }

    float
    Reduction<rt_max, tags::Cell>::value(const SparseVector<float> & a)
    {
        CONTEXT("When reducing SparseVector<float> to Scalar by maximum (Cell):");

        float result(0.0f);

        Operand oa = { &result };
        Operand ob = { a.elements() };
        Operand oc, od, oe;
        oc.u = a.used_elements() / 4096;
        od.u = a.used_elements() % 4096;
        od.u &= ~0xF;
        oe.f = a.used_elements() < a.size() ? 0.0f : a[*a.indices()];
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

        SPEInstruction instruction(oc_reduction_max_dense_float, 16 * 1024, oa, ob, oc, od, oe);

        if (use_spe)
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(oe.f);

        for (unsigned i(rest_index); i < a.used_elements() ; i++)
        {
            ppu_result = *(a.elements() + i) > ppu_result ? *(a.elements() + i) : ppu_result;
        }

        if (use_spe)
        {
            instruction.wait();
            result = ppu_result > result ? ppu_result : result;
        }
        else
        {
            result = ppu_result;
        }

        return result;
    }

    DenseVector<float>
    Reduction<rt_max, tags::Cell>::value(const SparseMatrix<float> & a)
    {
        CONTEXT("When reducing SparseMatrix<float> to Vector by max (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (SparseMatrix<float>::ConstRowIterator i(a.begin_non_zero_rows()), 
                i_end(a.end_non_zero_rows()) ; i != i_end ; ++i)
        {
            result[i.index()] = Reduction<rt_max, tags::Cell>::value(*i);
        }

        return result;
    }

}
