/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <cuda_util.hh>

extern "C" void * cuda_malloc(unsigned long bytes)
{
    void * gpu;
    cudaMalloc((void**)&gpu, bytes);
    CUDA_ERROR();
    return gpu;
}
extern "C" void cuda_upload(void * cpu, void * gpu, unsigned long bytes)
{
    cudaMemcpy(gpu, cpu, bytes, cudaMemcpyHostToDevice);
    CUDA_ERROR();
}

extern "C" void cuda_download(void * gpu, void * cpu, unsigned long bytes)
{
    cudaMemcpy(cpu, gpu, bytes, cudaMemcpyDeviceToHost);
    CUDA_ERROR();
}

extern "C" void cuda_free(void * gpu)
{
    cudaFree(gpu);
    CUDA_ERROR();
}

extern "C" void cuda_copy(void * src, void * dest, unsigned long bytes)
{
    cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToDevice);
    CUDA_ERROR();
}


