/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007,2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Till Barz <till.barz@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/cell/cell.hh>
#include <honei/la/reduction.hh>
#include <honei/util/memory_backend_cell.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_manager.hh>
#include <honei/util/stringify.hh>
#include <honei/util/partitioner.hh>

namespace honei
{
    using namespace cell;

    float
    Reduction<rt_sum, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<float> to Scalar by sum (Cell):");

        float ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::reduction_sum_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, float, rtm_mail> * > instructions;
        std::list<float> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            ppu_result += a[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                spu_results.push_back(0.0f);
                SPEFrameworkInstruction<1, float, rtm_mail> * instruction = new SPEFrameworkInstruction<1, float, rtm_mail>(
                        oc_reduction_sum_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, float, rtm_mail> * instruction = new SPEFrameworkInstruction<1, float, rtm_mail>(
                    oc_reduction_sum_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                ppu_result += a[index];
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<1, float, rtm_mail> * >::iterator i(instructions.begin()),
                    i_end(instructions.end()) ; i != i_end ; ++i)
            {
                if ((*i)->use_spe())
                    (*i)->wait();

                delete *i;
            }
        }
        for (std::list<float>::iterator i(spu_results.begin()), i_end(spu_results.end()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        return ppu_result;
    }

    DenseVector<float>
    Reduction<rt_sum, tags::Cell>::value(const DenseMatrix<float> & a)
    {
        /// \todo Add multi spu version after one spu programm can produce the whole vector.
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

    // double
    double
    Reduction<rt_sum, tags::Cell>::value(const DenseVectorContinuousBase<double> & a)
    {
        CONTEXT("When reducing DenseVectorContinuousBase<double> to Scalar by sum (Cell):");

        double ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x1);

        unsigned long spe_count(Configuration::instance()->get_value("cell::reduction_sum_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, double, rtm_mail> * > instructions;
        std::list<double> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            ppu_result += a[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                spu_results.push_back(0.0f);
                SPEFrameworkInstruction<1, double, rtm_mail> * instruction = new SPEFrameworkInstruction<1, double, rtm_mail>(
                        oc_reduction_sum_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, double, rtm_mail> * instruction = new SPEFrameworkInstruction<1, double, rtm_mail>(
                    oc_reduction_sum_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                ppu_result += a[index];
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<1, double, rtm_mail> * >::iterator i(instructions.begin()),
                    i_end(instructions.end()) ; i != i_end ; ++i)
            {
                if ((*i)->use_spe())
                    (*i)->wait();

                delete *i;
            }
        }
        for (std::list<double>::iterator i(spu_results.begin()), i_end(spu_results.end()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        return ppu_result;
    }

    DenseVector<double>
    Reduction<rt_sum, tags::Cell>::value(const DenseMatrix<double> & a)
    {
        /// \todo Add multi spu version after one spu programm can produce the whole vector.
        CONTEXT("When reducing DenseMatrix<double> to Vector by sum (Cell):");

        DenseVector<double> result(a.rows(), 0.0);

        for (unsigned k(0) ; k < a.rows() ; k++)
        {
            DenseVectorRange<double> row(a[k]);

            SPEFrameworkInstruction<1, double, rtm_mail> instruction(oc_reduction_sum_dense_double, result.elements() + k,
                    row.elements(), row.size());

            if (instruction.use_spe())
            {
                SPEManager::instance()->dispatch(instruction);
            }

            double ppu_result(0.0);

            for (Vector<double>::ConstElementIterator i(row.begin_elements()), i_end(row.element_at(instruction.transfer_begin())) ;
                    i != i_end ; ++i)
            {
                ppu_result += *i;
            }

            for (Vector<double>::ConstElementIterator i(row.element_at(instruction.transfer_end())), i_end(row.end_elements()) ; i != i_end ; ++i)
            {
                ppu_result += *i;
            }

            if (instruction.use_spe())
                instruction.wait();

            result[k] += ppu_result;

        }

        return result;
    }

    //

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

    //
    double
    Reduction<rt_sum, tags::Cell>::value(const SparseVector<double> & a)
    {
        CONTEXT("When reducing SparseVector<double> to Scalar by sum (Cell):");

        double result(0.0);
        SPEFrameworkInstruction<1, double, rtm_mail> instruction(oc_reduction_sum_dense_double, &result, a.elements(), a.used_elements());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        double ppu_result(0.0);

        for (Vector<double>::ConstElementIterator i(a.non_zero_element_at(instruction.transfer_end())),
                i_end(a.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            ppu_result += *i;
        }

        if (instruction.use_spe())
            instruction.wait();

        return result + ppu_result;
    }

    DenseVector<double>
    Reduction<rt_sum, tags::Cell>::value(const SparseMatrix<double> & a)
    {
        CONTEXT("When reducing SparseMatrix<double> to Vector by sum (Cell):");

        DenseVector<double> result(a.rows(), 0.0);

        for (SparseMatrix<double>::ConstRowIterator i(a.begin_non_zero_rows()),
                i_end(a.end_non_zero_rows()) ; i != i_end ; ++i)
        {
            result[i.index()] = Reduction<rt_sum, tags::Cell>::value(*i);
        }

        return result;
    }

    //
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

