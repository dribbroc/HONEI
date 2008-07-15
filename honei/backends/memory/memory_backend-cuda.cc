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

#include <honei/util/exception.hh>
#include <honei/backends/memory/memory_backend.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/backends/cuda/transfer.hh>
#include <honei/util/memory_arbiter.hh>

/// \todo How can we handle multiple device memories (numa wise)? -> One must know the choosen device (multiple memory_values).
namespace honei
{
    void * MemoryBackend<tags::GPU::CUDA>::upload(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When uploading data (CUDA):");
        std::map<void *, void *>::iterator i(_address_map.find(address));
        if (i == _address_map.end())
        {
            void * device(this->alloc(memid, address, bytes));
            cuda_upload(address, device, bytes);
            return device;
        }
        else
        {
            return i->second;
        }
    }

    void MemoryBackend<tags::GPU::CUDA>::download(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When downloading data (CUDA):");
        std::map<void *, void *>::iterator i(_address_map.find(address));
        if (i == _address_map.end())
        {
            throw InternalError("MemoryBackend<tags::GPU::CUDA> download address not found!");
        }
        else
        {
            cuda_download(i->second, address, bytes);
        }

    }

    void * MemoryBackend<tags::GPU::CUDA>::alloc(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When allocating data (CUDA):");
        void * device(cuda_malloc(bytes));
        _id_map.insert(std::pair<void *, Chunk>(memid, Chunk(address, device, bytes)));
        _address_map.insert(std::pair<void *, void *>(address, device));
        return device;
    }

    void MemoryBackend<tags::GPU::CUDA>::free(void * memid)
    {
        CONTEXT("When freeing data (CUDA):");
        std::multimap<void *, Chunk>::iterator i;
        std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
        for (std::multimap<void *, Chunk>::iterator i(range.first) ; i != range.second ; ++i)
        {
            cuda_free(i->second.device);
            _address_map.erase(i->second.address);
        }
        _id_map.erase(range.first, range.second);
    }

    template class InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>;
    template class MemoryBackend<tags::GPU::CUDA>;
}

using namespace honei;
static MemoryBackendRegistrator gpu_backend_registrator(tags::GPU::CUDA::memory_value, MemoryBackend<tags::GPU::CUDA>::instance());
