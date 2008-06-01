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

#ifndef CUDA_GUARD_CUDA_UTIL_HH
#define CUDA_GUARD_CUDA_UTIL_HH 1

#include <honei/util/assertion.hh>

#if defined (DEBUG)
namespace honei
{
    namespace cuda
    {
        static cudaError honeiCudaError;
    }
}
#endif

namespace honei
{

/**
 * \def ASSERT
 *
 * \brief Convenience definition that provides a way to test for cuda errors.
 *
 * The thrown Assertion will be automatically provided with the correct filename,
 * line number and function name.
 *
 * \warning Will only be compiled in when debug support is enabled.
 *
 * \ingroup grpcuda
 */
#if defined (DEBUG)
#define CUDA_ERROR() \
    do { \
        honei::cuda::honeiCudaError = cudaGetLastError(); \
        EXTERNAL_ASSERT(honei::cuda::honeiCudaError == cudaSuccess, cudaGetErrorString(honei::cuda::honeiCudaError)); \
    } while (false)
#else
#define CUDA_ERROR()
#endif
}

#endif
