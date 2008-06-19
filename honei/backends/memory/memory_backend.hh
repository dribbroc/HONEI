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

#ifndef MEMORY_GUARD_MEMORY_BACKEND_HH
#define MEMORY_GUARD_MEMORY_BACKEND_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/memory_backend_base.hh>
#include <honei/util/tags.hh>
#include <honei/util/exception.hh>

#include <map>

namespace honei
{
    template<typename Tag_>
    class MemoryBackend :
        public MemoryBackendBase
    {
    };

    template<>
    class MemoryBackend<tags::CPU> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>
    {
    public:
        virtual void * upload(unsigned long memid, void * address, unsigned long size)
        {
            CONTEXT("When uploading data:");
            return address;
        }

        virtual void download(unsigned long memid, void * address, unsigned long size)
        {
            CONTEXT("When downloading data:");
        }

        virtual void free(unsigned long memid, void * address, unsigned long size)
        {
        }

    };

    template<>
    class MemoryBackend<tags::GPU::CUDA> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>
    {
        private:
            std::map<void *, void *> address_map;

        public:
            virtual void * upload(unsigned long memid, void * address, unsigned long size);

            virtual void download(unsigned long memid, void * address, unsigned long size);

            virtual void free(unsigned long memid, void * address, unsigned long size);
    };
}
#endif
