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

extern "C" unsigned long cuda_upload(unsigned long src, unsigned long bytes)
{
    unsigned long gpu;
    cudaMalloc((void**)&gpu, bytes);
    cudaMemcpy((void*)&gpu, (void*)&src, bytes, cudaMemcpyHostToDevice);
    return gpu;
}

extern "C" void cuda_download(unsigned long src, unsigned long target, unsigned long bytes)
{
    cudaMemcpy((void*)&target, (void*)&src, bytes, cudaMemcpyDeviceToHost);
}

extern "C" void cuda_free(unsigned long src)
{
    cudaFree((void*)&src);
}

