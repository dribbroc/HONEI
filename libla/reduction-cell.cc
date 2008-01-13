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

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_sum_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return result += ppu_result;
    }

    DenseVector<float>
    Reduction<rt_sum, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by sum (Cell):");

        DenseVector<float> result(a.rows(), 0.0f);

        for (unsigned k(0) ; k < a.rows() ; k++)
        {
            DenseVectorRange<float> row(a[k]);

            SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_sum_dense_float, result.elements() + k,
                    row.elements(), row.size());

            if (instruction.use_spe())
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(0.0f);

            for (Vector<float>::ConstElementIterator i(row.begin_elements()), i_end(row.element_at(instruction.transfer_begin())) ;
                    i != i_end ; ++i)
            {
                ppu_result += *i;
            }

            for (Vector<float>::ConstElementIterator i(row.element_at(instruction.transfer_end())), i_end(row.end_elements()) ; i != i_end ; ++i)
            {
                ppu_result += *i;
            }

            if (instruction.use_spe())
                instruction.wait();

            result[k] += ppu_result;

        }

        return result;
    }

    float
    Reduction<rt_min, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by minimum (Cell):");

        float result(std::numeric_limits<float>::max());

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_min_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(std::numeric_limits<float>::max());

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result < *i) ? ppu_result : *i;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result < *i) ? ppu_result : *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return (result < ppu_result) ? result : ppu_result;
    }

    DenseVector<float>
    Reduction<rt_min, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by minimum (Cell):");

        DenseVector<float> result(a.rows(), std::numeric_limits<float>::max());

        for (unsigned k(0) ; k < a.rows() ; k++)
        {
            DenseVectorRange<float> row(a[k]);

            SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_min_dense_float, result.elements() + k,
                    row.elements(), row.size());

            if (instruction.use_spe())
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(std::numeric_limits<float>::max());

            for (Vector<float>::ConstElementIterator i(row.begin_elements()), i_end(row.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
            {
                ppu_result = (ppu_result < *i) ? ppu_result : *i;
            }

            for (Vector<float>::ConstElementIterator i(row.element_at(instruction.transfer_end())), i_end(row.end_elements()) ; i != i_end ; ++i)
            {
                ppu_result = (ppu_result < *i) ? ppu_result : *i;
            }

            if (instruction.use_spe())
            {
                instruction.wait();
            }

            result[k] = ppu_result < result[k] ? ppu_result : result[k];
        }

        return result;
    }

    float
    Reduction<rt_max, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by maximum (Cell):");

        float result(std::numeric_limits<float>::min());

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_max_dense_float, &result, a.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(std::numeric_limits<float>::min());

        for (Vector<float>::ConstElementIterator i(a.begin_elements()), i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result > *i) ? ppu_result : *i;
        }

        for (Vector<float>::ConstElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result > *i) ? ppu_result : *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return (ppu_result > result) ? ppu_result : result;
    }

    DenseVector<float>
    Reduction<rt_max, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When reducing DenseMatrix<float> to Vector by sum (Cell):");

        DenseVector<float> result(a.rows(), std::numeric_limits<float>::min());

        for (unsigned k(0) ; k < a.rows() ; k++)
        {
            DenseVectorRange<float> row(a[k]);

            SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_max_dense_float, result.elements() + k,
                    row.elements(), row.size());

            if (instruction.use_spe())
            {
                SPEManager::instance()->dispatch(instruction);
            }

            float ppu_result(std::numeric_limits<float>::min());

            for (Vector<float>::ConstElementIterator i(row.begin_elements()), i_end(row.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
            {
                ppu_result = (ppu_result > *i) ? ppu_result : *i;
            }

            for (Vector<float>::ConstElementIterator i(row.element_at(instruction.transfer_end())), i_end(row.end_elements()) ; i != i_end ; ++i)
            {
                ppu_result = (ppu_result > *i) ? ppu_result : *i;
            }

            if (instruction.use_spe())
            {
                instruction.wait();
                result[k] = ppu_result > result[k] ? ppu_result : result[k];
            }
            else
            {
               result[k] = ppu_result;
            }
        }

        return result;
    }

    float
    Reduction<rt_sum, tags::Cell>::value(const SparseVector<float> & a)
    {
        CONTEXT("When reducing SparseVector<float> to Scalar by sum (Cell):");

        float result(0.0f);
        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_sum_dense_float, &result, a.elements(), a.used_elements());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.non_zero_element_at(instruction.transfer_end())),
                i_end(a.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return result + ppu_result;
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

        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_min_dense_float, &result, a.elements(), a.used_elements());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.non_zero_element_at(instruction.transfer_end())),
                i_end(a.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result < *i) ? ppu_result : *i;
        }

        if (instruction.use_spe())
        {
            instruction.wait();
        }

        return (ppu_result < result) ? ppu_result : result;
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
        SPEFrameworkInstruction<1, float, rtm_mail> instruction(oc_reduction_max_dense_float, &result, a.elements(), a.used_elements());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        float ppu_result(0.0f);

        for (Vector<float>::ConstElementIterator i(a.non_zero_element_at(instruction.transfer_end())),
                i_end(a.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            ppu_result = (ppu_result > *i) ? ppu_result : *i;
        }

        if (instruction.use_spe())
        {
            instruction.wait();
        }

        return (ppu_result > result) ? ppu_result : result;
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

