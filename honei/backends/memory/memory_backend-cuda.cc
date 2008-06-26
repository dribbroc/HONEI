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
#include <honei/backends/cuda/transfer.hh>

namespace honei
{
    void * MemoryBackend<tags::GPU::CUDA>::upload(unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When uploading data:");
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

    void MemoryBackend<tags::GPU::CUDA>::download(unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When downloading data:");
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

    void * MemoryBackend<tags::GPU::CUDA>::alloc(unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When allocating data:");
        void * device(cuda_malloc(bytes));
        _id_map.insert(std::pair<unsigned long, Chunk>(memid, Chunk(address, device, bytes)));
        _address_map.insert(std::pair<void *, void *>(address, device));
        return device;
    }

    void MemoryBackend<tags::GPU::CUDA>::free(unsigned long memid)
    {
        CONTEXT("When freeing data:");
        std::multimap<unsigned long, Chunk>::iterator i;
        std::pair<std::multimap<unsigned long, Chunk>::iterator, std::multimap<unsigned long, Chunk>::iterator> range(_id_map.equal_range(memid));
        for (std::multimap<unsigned long, Chunk>::iterator i(range.first) ; i != range.second ; ++i)
        {
            cuda_free(i->second.device);
            _address_map.erase(i->second.address);
        }
        _id_map.erase(range.first, range.second);
    }
}
