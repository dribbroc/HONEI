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

#include <honei/util/unittest.hh>
#include <honei/util/memory_backend.hh>

using namespace honei;
using namespace tests;

template <typename Tag_>
class MemoryBackendQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        MemoryBackendQuickTest() :
            QuickTaggedTest<Tag_>("memory_backend_test")
        {
        }

        virtual void run() const
        {
            int data_array [10];
            for (int i(0) ; i < 10 ; ++i)
            {
                data_array[i] = i;
            }
            void * data = data_array;

            void * device(MemoryBackend<Tag_>::instance()->upload(0, data, 5 * sizeof(int)));
            if (device != data)
            {
                data_array[3] = -50;
            }
            MemoryBackend<Tag_>::instance()->download(0, data, 5 * sizeof(int));
            MemoryBackend<Tag_>::instance()->free(0);
            TEST_CHECK_EQUAL(data_array[3], 3);
        }
};
MemoryBackendQuickTest<tags::CPU> memory_backend_quick_test;
