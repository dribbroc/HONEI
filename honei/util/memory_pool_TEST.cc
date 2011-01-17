/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/unittest.hh>
#include <honei/util/memory_pool.hh>

using namespace honei;
using namespace tests;

template <typename Tag_>
class MemoryPoolQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        MemoryPoolQuickTest() :
            QuickTaggedTest<Tag_>("memory_pool_test")
        {
        }

        virtual void run() const
        {
            void * data = MemoryPool<Tag_>::instance()->alloc(10 * sizeof(int));
            int * data_array((int*)data);
            data_array[3] = 3;
            MemoryPool<Tag_>::instance()->free(data);
            data = MemoryPool<Tag_>::instance()->alloc(10 * sizeof(int));
            data_array = (int*)data;
            TEST_CHECK_EQUAL(data_array[3], 3);
            MemoryPool<Tag_>::instance()->free(data);
            MemoryPool<Tag_>::instance()->release_free();

            void * data2 = MemoryPool<Tag_>::instance()->alloc(30 * sizeof(unsigned long));
            data2 = MemoryPool<Tag_>::instance()->realloc(data2, 35 * sizeof(unsigned long));
            MemoryPool<Tag_>::instance()->free(data2);

        }
};
MemoryPoolQuickTest<tags::CPU> memory_pool_quick_test;
