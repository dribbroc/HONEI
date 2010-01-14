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

#include <honei/backends/cuda/cuda_util.hh>
#include <stdio.h>

extern "C" int cuda_device_count()
{
    int device_count;
    cudaGetDeviceCount(&device_count);
    CUDA_ERROR();
    return device_count;
}

extern "C" void cuda_print_device_name(int device)
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);
    printf("%s", prop.name);
}

extern "C" int cuda_get_device()
{
    int device(4711);
    cudaGetDevice(&device);
    return device;
}

extern "C" void cuda_set_device(int device)
{
    cudaSetDevice(device);
}
