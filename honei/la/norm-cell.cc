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

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/ppe/memory_backend_cell.hh>
#include <honei/backends/cell/ppe/spe_instruction.hh>
#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/la/norm.hh>
#include <honei/util/stringify.hh>
#include <honei/util/partitioner.hh>

#include <cmath>

namespace honei
{
    using namespace cell;

    float
    Norm<vnt_max, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        /// \todo Refactor to multi spu variant.
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

        float ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::norm_l_one_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, float, rtm_mail> * > instructions;
        std::list<float> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            ppu_result += fabs(a[index]);
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
                        oc_norm_l_one_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, float, rtm_mail> * instruction = new SPEFrameworkInstruction<1, float, rtm_mail>(
                    oc_norm_l_one_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                ppu_result += fabs(a[index]);
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

    double
    Norm<vnt_l_one, false, tags::Cell>::value(const DenseVectorContinuousBase<double> & a)
    {
        CONTEXT("When applying L1-norm to DenseVectorContinuousBase<double> (Cell):");

        double ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x1);

        unsigned long spe_count(Configuration::instance()->get_value("cell::norm_l_one_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, double, rtm_mail> * > instructions;
        std::list<double> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            ppu_result += fabs(a[index]);
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
                        oc_norm_l_one_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, double, rtm_mail> * instruction = new SPEFrameworkInstruction<1, double, rtm_mail>(
                    oc_norm_l_one_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                ppu_result += fabs(a[index]);
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

    float
    Norm<vnt_l_two, true, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm (with square-root) to DenseVectorContinuousBase<float> (Cell):");

        return sqrt(Norm<vnt_l_two, false, tags::Cell>::value(a));
    }

    double
    Norm<vnt_l_two, true, tags::Cell>::value(const DenseVectorContinuousBase<double> & a)
    {
        CONTEXT("When applying L2-norm (with square-root) to DenseVectorContinuousBase<double> (Cell):");

        return sqrt(Norm<vnt_l_two, false, tags::Cell>::value(a));
    }

    float
    Norm<vnt_l_two, false, tags::Cell>::value(const DenseVectorContinuousBase<float> & a)
    {
        CONTEXT("When applying L2-norm to DenseVectorContinuousBase<float> (Cell):");

        float ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::norm_l_two_dense_float", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, float, rtm_mail> * > instructions;
        std::list<float> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            float temp(fabs(a[index]));
            temp *= temp;
            ppu_result += temp;
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
                        oc_norm_l_two_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, float, rtm_mail> * instruction = new SPEFrameworkInstruction<1, float, rtm_mail>(
                    oc_norm_l_two_dense_float, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                float temp(fabs(a[index]));
                temp *= temp;
                ppu_result += temp;
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

    double
    Norm<vnt_l_two, false, tags::Cell>::value(const DenseVectorContinuousBase<double> & a)
    {
        CONTEXT("When applying L2-norm to DenseVectorContinuousBase<double> (Cell):");

        double ppu_result(0.0f);
        unsigned long skip(a.offset() & 0x1);

        unsigned long spe_count(Configuration::instance()->get_value("cell::norm_l_two_dense_double", 2ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<1, double, rtm_mail> * > instructions;
        std::list<double> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            double temp(fabs(a[index]));
            temp *= temp;
            ppu_result += temp;
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
                        oc_norm_l_two_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<1, double, rtm_mail> * instruction = new SPEFrameworkInstruction<1, double, rtm_mail>(
                    oc_norm_l_two_dense_double, &(spu_results.back()), a.elements() + skip + p->start, p->size);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                double temp(fabs(a[index]));
                temp *= temp;
                ppu_result += temp;
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
}

