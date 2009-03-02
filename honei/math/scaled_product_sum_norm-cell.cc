/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Danny van Dyk <dyk@honei.org>
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
#include <honei/backends/cell/ppe/spe_instruction.hh>
#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/math/scaled_product_sum_norm.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>

namespace honei
{
    using namespace cell;

    float
    ScaledProductSumNorm_TUTORIAL<tags::Cell>::value(float alpha, DenseVector<float> & y, float beta, BandedMatrixQ1<float> & A, DenseVector<float> & x)
    {
        CONTEXT("When calculating weighted norm of DenseVector<float> and DenseVector<float> (Cell):");


        //Still use HONEIs BandedMatrix-DenseVector product:
        DenseVector<float> b(Product<tags::Cell>::value(A, x));
        DenseVector<float> a(y);

        float ppu_result(0.0f);

        unsigned long skip(a.offset() & 0x3);
        if (0 != skip)
            skip = 4 - skip;

        unsigned long spe_count(Configuration::instance()->get_value("cell::scaled_product_sum_norm_float", 4ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());

        std::list<SPEFrameworkInstruction<2, float, rtm_mail> * > instructions;
        std::list<float> spu_results;
        PartitionList partitions;

        // Calculate the first elements on PPU (if needed).
        for (unsigned long index(0) ; index < skip && index < a.size() ; ++index)
        {
            float sum = alpha * a[index] + beta * b[index];
            ppu_result += sum * sum;
        }

        if (skip < a.size())
        {
            Partitioner<tags::Cell>(spe_count, std::max(a.size() / spe_count, 16ul), a.size() - skip, PartitionList::Filler(partitions));
            // Assemble instructions.
            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ;
                    p != p_last ; ++p)
            {
                spu_results.push_back(0.0f);
                SPEFrameworkInstruction<2, float, rtm_mail> * instruction = new SPEFrameworkInstruction<2, float, rtm_mail>(
                        oc_scaled_product_sum_norm_float, &(spu_results.back()), a.elements() + skip + p->start,
                        b.elements() + skip + p->start, p->size, alpha, beta);

                if (instruction->use_spe())
                {
                    SPEManager::instance()->dispatch(*instruction);
                    instructions.push_back(instruction);
                }
            }


            PartitionList::ConstIterator p(partitions.last());
            spu_results.push_back(0.0f);
            SPEFrameworkInstruction<2, float, rtm_mail> * instruction = new SPEFrameworkInstruction<2, float, rtm_mail>(
                    oc_scaled_product_sum_norm_float, &(spu_results.back()), a.elements() + skip + p->start,
                    b.elements() + skip + p->start, p->size, alpha, beta);

            if (instruction->use_spe())
            {
                SPEManager::instance()->dispatch(*instruction);
                instructions.push_back(instruction);
            }


            // Calculate the last elements on PPU (if needed).
            for (unsigned long index(skip + p->start + instruction->transfer_end()) ; index < a.size() ; ++index)
            {
                float sum = alpha * a[index] + beta * b[index];
                ppu_result += sum * sum;
            }

            // Wait for the SPU side
            for (std::list<SPEFrameworkInstruction<2, float, rtm_mail> * >::iterator i(instructions.begin()),
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

        return /*std::sqrt*/(ppu_result);
    }
}
