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

#include <iostream>
#include <honei/util/unittest.hh>
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/util/memory_backend.hh>

using namespace honei;
using namespace tests;

template <typename Tag_>
class MultiGPUQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        MultiGPUQuickTest() :
            QuickTaggedTest<Tag_>("multi_gpu_test")
        {
        }

        virtual void run() const
        {
            int device_count(-1);
            device_count = cuda_device_count();
            std::cout<<"Device Count: " << device_count<<std::endl;
            TEST_CHECK(device_count > 0);
            for (int device(0) ; device < device_count; ++device)
            {
                std::cout<<"Device " << device << ": " << std::endl;
                cuda_print_device_info(device);
            }

        }
};
MultiGPUQuickTest<tags::GPU::CUDA> cuda_multi_gpu_quick_test;
