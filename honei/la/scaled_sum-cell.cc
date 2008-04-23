/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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

#include <honei/backends/cell/cell.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/util/memory_backend_cell.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_manager.hh>
#include <honei/util/stringify.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/configuration.hh>

namespace honei
{
    using namespace cell;

    DenseVectorContinuousBase<float> &
    ScaledSum<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, const float & c)
    {
        CONTEXT("When adding DenseVectorContinuousBase<float> to scaled DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::scaled_sum_dense_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b[index] * c;
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                        oc_scaled_sum_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size, c);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                    oc_scaled_sum_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size, c);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b[index] * c;
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<2, float, rtm_dma> * >::iterator i(instructions.begin()),
                    i_end(instructions.end()) ; i != i_end ; ++i)
            {
                if ((*i)->use_spe())
                    (*i)->wait();

                delete *i;
            }
        }

        return a;
    }

    DenseVectorContinuousBase<double> &
    ScaledSum<tags::Cell>::value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, const double & c)
    {
        CONTEXT("When adding DenseVectorContinuousBase<double> to scaled DenseVectorContinuousBase<double> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x1);

        unsigned long spe_count(Configuration::instance()->get_value("cell::scaled_sum_dense_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b[index] * c;
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                        oc_scaled_sum_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size, c);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                    oc_scaled_sum_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size, c);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b[index] * c;
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<2, double, rtm_dma> * >::iterator i(instructions.begin()),
                    i_end(instructions.end()) ; i != i_end ; ++i)
            {
                if ((*i)->use_spe())
                    (*i)->wait();

                delete *i;
            }
        }

        return a;
    }

    DenseMatrix<float> &
    ScaledSum<tags::Cell>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b, const float & c)
    {
        CONTEXT("When adding DenseMatrix<float> to scaled DenseMatrix<float> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        SPEFrameworkInstruction<2, float, rtm_dma> instruction(oc_scaled_sum_dense_dense_float, a.elements(), b.elements(), a.rows() * a.columns(), c);

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        // Calculate the first 4 - a_offset % 4 elements on PPU.
        Matrix<float>::ConstElementIterator j(b.begin_elements());
        for (MutableMatrix<float>::ElementIterator i(a.begin_elements()),
            i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i, ++j)
        {
            *i += *j * c;
        }

        Matrix<float>::ConstElementIterator k(b.element_at(instruction.transfer_end()));
        for (MutableMatrix<float>::ElementIterator i(a.element_at(instruction.transfer_end())),
            i_end(a.end_elements()) ; i != i_end ; ++i, ++k)
        {
            *i += *k * c;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }
}
