/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <unittest/unittest.hh>
#include <honei/util/memory_arbiter.hh>

using namespace honei;
using namespace tests;

class MemoryArbiterQuickTest :
    public QuickTest
{
    public:
        MemoryArbiterQuickTest() :
            QuickTest("memory_arbiter_test")
        {
        }

        virtual void run() const
        {
            char  data_array1 [10];
            void * data1 = data_array1;
            void * mem1(data1);
            char  data_array2 [10];
            void * data2 = data_array2;
            void * mem2(data2);

            MemoryArbiter::instance()->register_address(mem2);
            MemoryArbiter::instance()->register_address(mem1);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem1, data1, 2);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->unlock(lm_read_only, mem2);
            MemoryArbiter::instance()->unlock(lm_read_only, mem1);
            MemoryArbiter::instance()->lock(lm_read_and_write, tags::CPU::memory_value, mem1, data1, 2);
            MemoryArbiter::instance()->unlock(lm_read_and_write, mem1);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->unlock(lm_read_only, (void *)25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->unlock(lm_read_only, mem1), InternalError);

            MemoryArbiter::instance()->remove_address(mem2);
            MemoryArbiter::instance()->remove_address(mem1);
        }
} memory_arbiter_quick_test;
