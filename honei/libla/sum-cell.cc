/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/libla/sum.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/partitioner.hh>

#include <iostream>

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

        unsigned long spe_count(std::min(2ul, SPEManager::instance()->spe_count()));

        std::list<SPEFrameworkInstruction<2, float, rtm_dma> * > instructions;
        PartitionList partitions;
        Partitioner(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));

        // Assemble instructions.
        for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                p != p_last ; ++p)
        {
            SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(oc_sum_dense_dense_float,
                    a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
            }

            instructions.push_back(instruction);

        }

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip ; ++index)
        {
            a[index] += b[index];
        }

        // Calculate the last elements on PPU (if needed).
        PartitionList::ConstIterator p(partitions.last());
        for (unsigned long index(skip + p->start) ; index < skip + p->start + p->size ; ++index)
        {
            a[index] += b[index];
        }

        // Wait for the SPU side
        for (std::list<SPEFrameworkInstruction<2, float, rtm_dma> * >::iterator i(instructions.begin()), i_end(instructions.end()) ;
                i != i_end ; ++i)
        {
            if ((*i)->use_spe())
                (*i)->wait();
            delete *i;
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

        SPEFrameworkInstruction<2, float, rtm_dma> instruction(oc_sum_dense_dense_float, a.elements(), b.elements(), a.rows() * a.columns());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        Matrix<float>::ConstElementIterator k(b.element_at(instruction.transfer_end()));
        for (MutableMatrix<float>::ElementIterator i(a.element_at(instruction.transfer_end())),
            i_end(a.end_elements()) ; i != i_end ; ++i, ++k)
        {
            *i += *k;
        }

        if (instruction.use_spe())
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

        SPEFrameworkInstruction<1, float, rtm_dma> instruction(oc_sum_dense_scalar_float, a.elements(), a.size(), b);

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (Vector<float>::ElementIterator i(a.begin_elements()),
            i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i)
        {
            *i += b;
        }

        for (Vector<float>::ElementIterator i(a.element_at(instruction.transfer_end())),
                i_end(a.end_elements()) ; i != i_end ; ++i)
        {
            *i += b;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }

    DenseMatrix<float> &
    Sum<tags::Cell>::value(DenseMatrix<float> & a, const float & b)
    {
        CONTEXT("When adding float to DenseMatrix<float> (Cell):");

        SPEFrameworkInstruction<1, float, rtm_dma> instruction(oc_sum_dense_scalar_float, a.elements(), a.rows() * a.columns(), b);

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        for (MutableMatrix<float>::ElementIterator i(a.element_at(instruction.transfer_end())), i_end(a.end_elements()) ;
                i != i_end ; ++i)
        {
            *i += b;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }
}
