/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
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
    unsigned long MemoryBackend<tags::GPU::CUDA>::upload(unsigned long memid, unsigned long address, unsigned long bytes)
    {
        std::map<unsigned long, unsigned long>::iterator i(address_map.find(address));
        if (i == address_map.end())
        {
            unsigned long temp(cuda_upload(address, bytes));
            address_map.insert(std::pair<unsigned long, unsigned long>(address, temp));
            return temp;
        }
        else
        {
            return i->second;
        }
    }

    void MemoryBackend<tags::GPU::CUDA>::download(unsigned long memid, unsigned long address, unsigned long bytes)
    {
        std::map<unsigned long, unsigned long>::iterator i(address_map.find(address));
        if (i == address_map.end())
        {
            throw InternalError("MemoryBackend<tags::GPU::CUDA> download address not found!");
        }
        else
        {
            cuda_download(i->second, address, bytes);
            cuda_free(i->second);
            address_map.erase(i);
        }

    }

    void MemoryBackend<tags::GPU::CUDA>::reset(unsigned long memid, unsigned long address, unsigned long size)
    {
        std::map<unsigned long, unsigned long>::iterator i(address_map.find(address));
        if (i != address_map.end())
        {
            cuda_free(i->second);
            address_map.erase(i);
        }
    }
}
