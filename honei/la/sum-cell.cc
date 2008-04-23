/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/la/sum.hh>
#include <honei/util/memory_backend_cell.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_manager.hh>
#include <honei/util/stringify.hh>
#include <honei/util/partitioner.hh>

namespace honei
{
    using namespace cell;

    DenseVectorContinuousBase<float> &
    Sum<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When adding DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                        oc_sum_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                    oc_sum_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b[index];
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
    Sum<tags::Cell>::value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b)
    {
        CONTEXT("When adding DenseVectorContinuousBase<double> to DenseVectorContinuousBase<double> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x1);
        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                        oc_sum_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                    oc_sum_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b[index];
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
    Sum<tags::Cell>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When adding DenseMatrix<float> to DenseMatrix<float> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).

        Partitioner<tags::Cell>(spe_count, std::max(a.rows() * a.columns() / spe_count, 16ul), a.rows() * a.columns(), PartitionList::Filler(partitions));
        // Assemble instructions.
        for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                p != p_last ; ++p)
        {
            SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                    oc_sum_dense_dense_float, a.elements() + p->start, b.elements() + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                oc_sum_dense_dense_float, a.elements() + p->start, b.elements() + p->start, p->size);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] += b.elements()[index];
        }

        // Wait for the SPU side
        for (std::list<SPEFrameworkInstruction<2, float, rtm_dma> * >::iterator i(instructions.begin()),
                i_end(instructions.end()) ; i != i_end ; ++i)
        {
            if ((*i)->use_spe())
                (*i)->wait();

            delete *i;
        }


        return a;
    }

    DenseMatrix<double> &
    Sum<tags::Cell>::value(DenseMatrix<double> & a, const DenseMatrix<double> & b)
    {
        CONTEXT("When adding DenseMatrix<double> to DenseMatrix<double> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).

        Partitioner<tags::Cell>(spe_count, std::max(a.rows() * a.columns() / spe_count, 16ul), a.rows() * a.columns(), PartitionList::Filler(partitions));
        // Assemble instructions.
        for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                p != p_last ; ++p)
        {
            SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                    oc_sum_dense_dense_double, a.elements() + p->start, b.elements() + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                oc_sum_dense_dense_double, a.elements() + p->start, b.elements() + p->start, p->size);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] += b.elements()[index];
        }

        // Wait for the SPU side
        for (std::list<SPEFrameworkInstruction<2, double, rtm_dma> * >::iterator i(instructions.begin()),
                i_end(instructions.end()) ; i != i_end ; ++i)
        {
            if ((*i)->use_spe())
                (*i)->wait();

            delete *i;
        }


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

        unsigned long last_index(b.indices()[b.used_elements()-1]);

        unsigned a_needed_transfers = ((last_index + 1) * 4) / 16384;
        unsigned a_rest_transfer_size(((last_index + 1) * 4) % 16384);
        if (a_rest_transfer_size != 0)
            a_needed_transfers++;

        if (a_needed_transfers == 0)
            a_needed_transfers++;

        // We can only transfer half of the B-part, cause indices pointer is 64 bit.
        unsigned long b_needed_transfers((od.u * 4) / 8192);
        unsigned b_rest_transfer_size((od.u * 4) % 8192);
        if (b_rest_transfer_size != 0)
            b_needed_transfers++;

        if (b_needed_transfers == 0)
            b_needed_transfers++;

        Operand oe, of, og, oh, oi;

        oe.u = a_needed_transfers;
        of.u = a_rest_transfer_size;
        og.u = b_needed_transfers;
        oh.u = b_rest_transfer_size;
        // Last index of first transfered elements of A
        oi.u = a.size() > 4096 ? 4095 : a.size()-1;

        Operand oj;
        oj.u = a.size();

        SPEInstruction instruction(oc_sum_dense_sparse_float, 16384, oa, ob, oc, od, oe, of, og, oh, oi, oj);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();

        return a;
    }

    DenseVectorContinuousBase<float> &
    Sum<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const float & b)
    {
        CONTEXT("When adding float to DenseVectorContinuousBase<float> (Cell):");

        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_scalar_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b;
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<1, float, rtm_dma> * instruction = new SPEFrameworkInstruction<1, float, rtm_dma>(
                        oc_sum_dense_scalar_float, a.elements() + skip + p->start, p->size, b);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<1, float, rtm_dma> * instruction = new SPEFrameworkInstruction<1, float, rtm_dma>(
                    oc_sum_dense_scalar_float, a.elements() + skip + p->start, p->size, b);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b;
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<1, float, rtm_dma> * >::iterator i(instructions.begin()),
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
    Sum<tags::Cell>::value(DenseVectorContinuousBase<double> & a, const double & b)
    {
        CONTEXT("When adding double to DenseVectorContinuousBase<double> (Cell):");

        unsigned long skip(a.offset() & 0x1);

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_scalar_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] += b;
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<1, double, rtm_dma> * instruction = new SPEFrameworkInstruction<1, double, rtm_dma>(
                        oc_sum_dense_scalar_double, a.elements() + skip + p->start, p->size, b);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<1, double, rtm_dma> * instruction = new SPEFrameworkInstruction<1, double, rtm_dma>(
                    oc_sum_dense_scalar_double, a.elements() + skip + p->start, p->size, b);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] += b;
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<1, double, rtm_dma> * >::iterator i(instructions.begin()),
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
    Sum<tags::Cell>::value(DenseMatrix<float> & a, const float & b)
    {
        CONTEXT("When adding float to DenseMatrix<float> (Cell):");

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_scalar_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).

        Partitioner<tags::Cell>(spe_count, std::max(a.rows() * a.columns() / spe_count, 16ul), a.rows() * a.columns(), PartitionList::Filler(partitions));
        // Assemble instructions.
        for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                p != p_last ; ++p)
        {
            SPEFrameworkInstruction<1, float, rtm_dma> * instruction = new SPEFrameworkInstruction<1, float, rtm_dma>(
                    oc_sum_dense_scalar_float, a.elements() + p->start, p->size, b);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<1, float, rtm_dma> * instruction = new SPEFrameworkInstruction<1, float, rtm_dma>(
                oc_sum_dense_scalar_float, a.elements() + p->start, p->size, b);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] += b;
        }

        // Wait for the SPU side
        for (std::list<SPEFrameworkInstruction<1, float, rtm_dma> * >::iterator i(instructions.begin()),
                i_end(instructions.end()) ; i != i_end ; ++i)
        {
            if ((*i)->use_spe())
                (*i)->wait();

            delete *i;
        }

        return a;
    }

    DenseMatrix<double> &
    Sum<tags::Cell>::value(DenseMatrix<double> & a, const double & b)
    {
        CONTEXT("When adding double to DenseMatrix<double> (Cell):");

        unsigned long spe_count(Configuration::instance()->get_value("cell::sum_dense_scalar_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).

        Partitioner<tags::Cell>(spe_count, std::max(a.rows() * a.columns() / spe_count, 16ul), a.rows() * a.columns(), PartitionList::Filler(partitions));
        // Assemble instructions.
        for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                p != p_last ; ++p)
        {
            SPEFrameworkInstruction<1, double, rtm_dma> * instruction = new SPEFrameworkInstruction<1, double, rtm_dma>(
                    oc_sum_dense_scalar_double, a.elements() + p->start, p->size, b);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<1, double, rtm_dma> * instruction = new SPEFrameworkInstruction<1, double, rtm_dma>(
                oc_sum_dense_scalar_double, a.elements() + p->start, p->size, b);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] += b;
        }

        // Wait for the SPU side
        for (std::list<SPEFrameworkInstruction<1, double, rtm_dma> * >::iterator i(instructions.begin()),
                i_end(instructions.end()) ; i != i_end ; ++i)
        {
            if ((*i)->use_spe())
                (*i)->wait();

            delete *i;
        }

        return a;
    }
}
