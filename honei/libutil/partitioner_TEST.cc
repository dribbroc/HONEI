/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2008 Joachim Messer <joachim.messer@uni-dortmund.de>
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/partitioner.hh>
#include <unittest/unittest.hh>

#include <list>

#include <iostream> // <<-- Remove!

using namespace honei;
using namespace tests;

template <typename Tag_, unsigned long best_part_size_, unsigned long quantisation_>
class PartitionerTestCell :
    public QuickTest
{
    public:
        PartitionerTestCell() :
            QuickTest("partitioner_test<" + stringify(Tag_::name) + "," + stringify(best_part_size_) + "," + stringify(quantisation_) + ">")
        {
        }

        virtual void run_test(unsigned max_count, unsigned long overall_size) const
        {
            PartitionList partitions;
            Partitioner<Tag_>(max_count, best_part_size_, overall_size, PartitionList::Filler(partitions));

            unsigned long count(partitions.size());
            unsigned long sum(0);

            for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ; p != p_last ; ++p)
            {
                unsigned long partition_size(p->size);
                sum += partition_size;

                TEST_CHECK(0 == partition_size % quantisation_);
                TEST_CHECK(0 == sum % quantisation_);

                if (overall_size > best_part_size_)
                {
                    TEST_CHECK(partition_size >= best_part_size_);
                }
            }

            sum += partitions.last()->size;

            if (count > max_count)
            {
                TEST_CHECK(partitions.last()->size < 16);
            }

            TEST_CHECK(count <= max_count + 1);
            TEST_CHECK_EQUAL(sum, overall_size);
        }

        virtual void run() const
        {
            for (unsigned long j(1), j_end(33) ; j != j_end ; ++j)
            {
                for (unsigned long k(1), k_end(best_part_size_) ; k < k_end ; k += 10)
                {
                    run_test(j, k);
                }

                for (unsigned long k(best_part_size_), k_end(64ul*64*64*128) ; k < k_end ; k += 30101)
                {
                    run_test(j, k);
                }
            }
        }
};

PartitionerTestCell<tags::Cell, 16384, 16> partitioner_test_cell_16k_16;


template <typename Tag_, unsigned long best_part_size_, unsigned long quantisation_>
class PartitionerTestMC :
    public QuickTest
{
    public:
        PartitionerTestMC() :
            QuickTest("partitioner_test<" + stringify(Tag_::name) + "," + stringify(best_part_size_) + "," + stringify(quantisation_) + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long j(1), j_end(33) ; j != j_end ; ++j)
            {
                unsigned long max_count(j);

                for (unsigned long k(best_part_size_ << 1), k_end(1 << 16) ; k < k_end ; k += 100)
                {
                    unsigned long overall_size(k);

                    PartitionList partitions;
                    Partitioner<Tag_>(max_count, best_part_size_, quantisation_, overall_size, PartitionList::Filler(partitions));

                    unsigned long count(partitions.size());
                    unsigned long sum(0);

                    for (PartitionList::ConstIterator p(partitions.begin()), p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        unsigned long partition_size(p->size);
                        sum += partition_size;

                        TEST_CHECK_EQUAL(0, partition_size % quantisation_);

                        if (overall_size > best_part_size_)
                        {
                            TEST_CHECK(partition_size >= best_part_size_);
                        }
                    }

                    sum += partitions.last()->size;

                    TEST_CHECK(count <= max_count + 1);
                    TEST_CHECK_EQUAL(overall_size, sum);
                }
            }
        }
};

PartitionerTestMC<tags::CPU::MultiCore, 1024, 16> partitioner_test_mc_1k_16;
PartitionerTestMC<tags::CPU::MultiCore, 4096, 16> partitioner_test_mc_4k_16;
PartitionerTestMC<tags::CPU::MultiCore, 131072, 16> partitioner_test_mc_128k_16;

