/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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
#include <honei/libla/element_product.hh>
#include <honei/util/memory_backend_cell.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_manager.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float> &
    ElementProduct<tags::Cell>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When multiplying DenseMatrix<float> to DenseMatrix<float> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        unsigned long spe_count(Configuration::instance()->get_value("cell::element_product_dense_dense_float", 2ul));
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
                    oc_element_product_dense_dense_float, a.elements() + p->start, b.elements() + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                oc_element_product_dense_dense_float, a.elements() + p->start, b.elements() + p->start, p->size);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] *= b.elements()[index];
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
    ElementProduct<tags::Cell>::value(DenseMatrix<double> & a, const DenseMatrix<double> & b)
    {
        CONTEXT("When multiplying DenseMatrix<double> to DenseMatrix<double> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        unsigned long spe_count(Configuration::instance()->get_value("cell::element_product_dense_dense_double", 2ul));
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
                    oc_element_product_dense_dense_double, a.elements() + p->start, b.elements() + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }
        }


        PartitionList::ConstIterator p(partitions.last());
        SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                oc_element_product_dense_dense_double, a.elements() + p->start, b.elements() + p->start, p->size);

        if (instruction->use_spe())
        {
            SPEManager::instance()->dispatch(*instruction);
            instructions.push_back(instruction);
        }


        // Calculate the last elements on PPU (if needed).
        for (unsigned long index(p->start + instruction->transfer_end()) ; index < a.rows() * a.columns() ; ++index)
        {
            a.elements()[index] *= b.elements()[index];
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

    DenseVectorContinuousBase<float> &
    ElementProduct<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When multiplying DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::element_product_dense_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, float, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] *= b[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                        oc_element_product_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, float, rtm_dma> * instruction = new SPEFrameworkInstruction<2, float, rtm_dma>(
                    oc_element_product_dense_dense_float, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] *= b[index];
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
    ElementProduct<tags::Cell>::value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b)
    {
        CONTEXT("When multiplying DenseVectorContinuousBase<double> to DenseVectorContinuousBase<double> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        unsigned long skip(a.offset() & 0x1);
        unsigned long spe_count(Configuration::instance()->get_value("cell::element_product_dense_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, double, rtm_dma> * > instructions;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            a[index] *= b[index];
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                        oc_element_product_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            SPEFrameworkInstruction<2, double, rtm_dma> * instruction = new SPEFrameworkInstruction<2, double, rtm_dma>(
                    oc_element_product_dense_dense_double, a.elements() + skip + p->start, b.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                a[index] *= b[index];
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

}
