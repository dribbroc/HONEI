/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Joachim Messer <joachim.messer@uni-dortmund.de>
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

using namespace honei;
using namespace tests;

typedef std::list<std::pair<unsigned long, unsigned long> > PartitionList;

struct PartitionerTester
{
    public: 
        PartitionList & list;

        PartitionerTester(PartitionList & l) :list(l)
        {
        }

        void operator() (unsigned long start, unsigned long size)
        {
            list.push_back(std::make_pair(start, size));
        }

};

class PartitionerTest :
    public QuickTest
{
    public:
        PartitionerTest() :
            QuickTest("partitioner_test")
        {
        }

        virtual void run() const
        {
            unsigned long best_part_size(16 * 1024);
            unsigned long max_count(16);

            PartitionList list;
            PartitionerTester tester(list);
            unsigned long offset, count, ps, temp;
            bool aligned;

            for (unsigned long size(100); size < (1 << 19); size += 100)
            {
                Partitioner(max_count, best_part_size, size, tester);

                ps = tester.list.front().second;
                aligned = !(ps % 16);
                TEST_CHECK_EQUAL(aligned, true);

                if (size >= best_part_size)
                {
                    TEST_CHECK(ps >= best_part_size);
                }

                offset = 0;
                count = 0;

                for ( ; offset <= size - ps; offset += ps)
                {
                    tester.list.pop_front();
                    ++count;
                }

                TEST_CHECK(count <= max_count);

                temp = count * ps;

                if (size > offset)
                {
                    temp += tester.list.front().second;
                    tester.list.pop_front();
                }

                TEST_CHECK_EQUAL(size, temp);
            }
        }
} partitioner_test;
