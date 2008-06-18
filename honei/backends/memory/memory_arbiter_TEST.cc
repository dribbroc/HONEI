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
#include <honei/backends/memory/memory_arbiter.hh>

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
            MemoryArbiter::instance()->read<tags::CPU>(0, 1, 2);
            MemoryArbiter::instance()->read<tags::CPU>(1, 2 , 1);
            MemoryArbiter::instance()->release_read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(0);
            MemoryArbiter::instance()->write<tags::CPU>(0, 7, 2);
            MemoryArbiter::instance()->release_write<tags::CPU>(0);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(0), InternalError);
        }
} memory_arbiter_quick_test;


#ifdef HONEI_CUDA
class CUDAMemoryArbiterQuickTest :
    public QuickTest
{
    public:
        CUDAMemoryArbiterQuickTest() :
            QuickTest("cuda_memory_arbiter_test")
        {
        }

        virtual void run() const
        {
            MemoryArbiter::instance()->read<tags::CPU>(0, 1, 1);
            MemoryArbiter::instance()->read<tags::GPU::CUDA>(0, 1, 1);
            MemoryArbiter::instance()->read<tags::CPU>(1, 1, 1);
            MemoryArbiter::instance()->release_read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(0);
            MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(0);
            MemoryArbiter::instance()->write<tags::CPU>(0, 1, 1);
            MemoryArbiter::instance()->release_write<tags::CPU>(0);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(0), InternalError);
        }
} cuda_memory_arbiter_quick_test;
#endif
