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

#include <honei/backends/memory/memory_backend.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/memory_arbiter.hh>

namespace honei
{
    template class InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>;
    template class MemoryBackend<tags::CPU>;

#ifdef HONEI_CUDA
    template class InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>;
    template class MemoryBackend<tags::GPU::CUDA>;
#endif
}

using namespace honei;
static MemoryBackendRegistrator cpu_backend_registrator(tags::CPU::memory_value, MemoryBackend<tags::CPU>::instance());
#ifdef HONEI_CUDA
static MemoryBackendRegistrator gpu_backend_registrator(tags::GPU::CUDA::memory_value, MemoryBackend<tags::GPU::CUDA>::instance());
#endif
