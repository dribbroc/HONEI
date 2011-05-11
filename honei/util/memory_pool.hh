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

#pragma once
#ifndef UTIL_GUARD_MEMORY_POOL_HH
#define UTIL_GUARD_MEMORY_POOL_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/tags.hh>
#include <honei/util/mutex.hh>
#include <honei/util/lock.hh>
#include <honei/util/exception.hh>

#include <cstring>
#include <map>
#include <string>
#include <cstdlib>
#include <new>
#include <iostream>

namespace honei
{
    namespace intern
    {
#ifndef DOXYGEN
        static inline unsigned multiple_of_sixteen(unsigned v)
        {
            unsigned r(v & 0xF);
            return (r > 0) ? v + 16 - r : v;
        }
#endif
    }

    template<typename Tag_>
        class MemoryPool
        {
        };

    template<>
        class MemoryPool<tags::CPU> :
        public InstantiationPolicy<MemoryPool<tags::CPU>, Singleton>
        {
        private:
            /// Our mutex.
            Mutex * const _mutex;

            std::multimap<unsigned long, void *> _free_chunks;
            std::map<void*, unsigned long> _used_chunks;

            MemoryPool() :
                _mutex(new Mutex)
            {
            }

            ~MemoryPool()
            {
                {
                    Lock l(*_mutex);
                    //std::cout<<"free chunks size: "<<_free_chunks.size()<<std::endl;
                    //std::cout<<"used chunks size: "<<_used_chunks.size()<<std::endl;
                    _release_free();
                }
                if (_used_chunks.size() != 0)
                    throw InternalError("MemoryPool destructor called with elements still unfreed");
                delete _mutex;
            }

            void _release_free()
            {
                std::multimap<unsigned long, void*>::iterator free_it;
                for (free_it = _free_chunks.begin() ; free_it != _free_chunks.end() ; ++free_it)
                {
                    ::free(free_it->second);
                }
                _free_chunks.clear();
            }

        public:
            friend class InstantiationPolicy<MemoryPool, Singleton>;

            void * alloc(unsigned long bytes)
            {
                CONTEXT("When allocating data (CPU):");
                unsigned long real_size(intern::multiple_of_sixteen(bytes));
                /*void * result(0);
                posix_memalign(&result, 16, real_size);
                return result;*/

                Lock l(*_mutex);
                std::multimap<unsigned long, void*>::iterator free_it;

                free_it = _free_chunks.find(real_size);
                if (free_it != _free_chunks.end())
                {
                    void * ret(free_it->second);
                    _used_chunks.insert (std::pair<void*, unsigned long>(free_it->second, free_it->first));
                    _free_chunks.erase(free_it);
                    return ret;
                }
                else
                {
                    void * result(0);
                    int status(0);

                    status = posix_memalign(&result, 16, real_size);
                    if (status == 0)
                    {
                        _used_chunks.insert (std::pair<void*, unsigned long>(result, real_size));
                        return result;
                    }
                    else
                    {
                        _release_free();
                        status = posix_memalign(&result, 16, real_size);
                    }
                    if (status == 0)
                    {
                        _used_chunks.insert (std::pair<void*, unsigned long>(result, real_size));
                        return result;
                    }
                    else
                        throw InternalError("MemoryPool: bad alloc or out of memory!");
                }

            }

            void * realloc(void * address, unsigned long bytes)
            {
                CONTEXT("When allocating data (CPU):");
                unsigned long real_size(intern::multiple_of_sixteen(bytes));

                Lock l(*_mutex);
                std::map<void *, unsigned long>::iterator used_it;

                used_it = _used_chunks.find(address);
                if (used_it != _used_chunks.end())
                {
                    if (real_size == used_it->second)
                    {
                        return address;
                    }
                    else
                    {
                        void * result(0);
                        int status(0);

                        status = posix_memalign(&result, 16, real_size);
                        if (status == 0)
                        {
                            _used_chunks.insert (std::pair<void*, unsigned long>(result, real_size));
                            //return result;
                        }
                        else
                        {
                            _release_free();
                            status = posix_memalign(&result, 16, real_size);
                        }
                        if (status == 0)
                        {
                            _used_chunks.insert (std::pair<void*, unsigned long>(result, real_size));
                            //return result;
                        }
                        else
                            throw InternalError("MemoryPool: bad alloc or out of memory!");

                        memcpy(result, address, used_it->second);
                        ::free(address);
                        _used_chunks.erase(used_it);
                        return result;
                    }
                }
                else
                {
                    throw InternalError("MemoryPool: realloc address not found!");
                }
            }

            void free(void * memid)
            {
                CONTEXT("When freeing data (CPU):");
                //::free(memid);
                Lock l(*_mutex);
                std::map<void*, unsigned long>::iterator used_it;

                used_it = _used_chunks.find(memid);
                if (used_it != _used_chunks.end())
                {
                    /// \todo reactivate free_chunk insertion, to store allocated chunks
                    //_free_chunks.insert(std::pair<unsigned long, void*>(used_it->second, used_it->first));
                    ::free(memid);
                    _used_chunks.erase(used_it);
                }
                else
                {
                    throw InternalError("MemoryPool: memory chunk not found!");
                }

            }

            void release_free()
            {
                CONTEXT("When releasing all memory chunks (CPU):");
                Lock l(*_mutex);
                _release_free();
            }
        };

    /*template<>
      class MemoryBackend<tags::GPU::CUDA> :
      public MemoryBackendBase,
      public InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>,
      public PrivateImplementationPattern<MemoryBackend<tags::GPU::CUDA>, Single>
      {
      public:
    /// Constructor.
    MemoryBackend<tags::GPU::CUDA>();

    /// Destructor.
    virtual ~MemoryBackend<tags::GPU::CUDA>();

    virtual void * upload(void * memid, void * address, unsigned long bytes);

    virtual void download(void * memid, void * address, unsigned long bytes);

    virtual void * alloc(void * memid, void * address, unsigned long bytes);

    virtual void free(void * memid);

    virtual void copy(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

    virtual bool knows(void * memid, void * address);
    };

    template<>
    class MemoryBackend<tags::OpenCL::CPU> :
    public MemoryBackendBase,
    public InstantiationPolicy<MemoryBackend<tags::OpenCL::CPU>, Singleton>,
    public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::CPU>, Single>
    {
    public:
    /// Constructor.
    MemoryBackend<tags::OpenCL::CPU>();

    /// Destructor.
    virtual ~MemoryBackend<tags::OpenCL::CPU>();

    virtual void * upload(void * memid, void * address, unsigned long bytes);

    virtual void download(void * memid, void * address, unsigned long bytes);

    virtual void * alloc(void * memid, void * address, unsigned long bytes);

    virtual void free(void * memid);

    virtual void copy(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
    void * dest_address, unsigned long bytes);

    virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

    virtual bool knows(void * memid, void * address);
    };

    template<>
    class MemoryBackend<tags::OpenCL::GPU> :
        public MemoryBackendBase,
               public InstantiationPolicy<MemoryBackend<tags::OpenCL::GPU>, Singleton>,
               public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::GPU>, Single>
               {
                   public:
                       /// Constructor.
                       MemoryBackend<tags::OpenCL::GPU>();

                       /// Destructor.
                       virtual ~MemoryBackend<tags::OpenCL::GPU>();

                       virtual void * upload(void * memid, void * address, unsigned long bytes);

                       virtual void download(void * memid, void * address, unsigned long bytes);

                       virtual void * alloc(void * memid, void * address, unsigned long bytes);

                       virtual void free(void * memid);

                       virtual void copy(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

                       virtual bool knows(void * memid, void * address);
               };

    template<>
        class MemoryBackend<tags::OpenCL::Accelerator> :
        public MemoryBackendBase,
               public InstantiationPolicy<MemoryBackend<tags::OpenCL::Accelerator>, Singleton>,
               public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::Accelerator>, Single>
               {
                   public:
                       /// Constructor.
                       MemoryBackend<tags::OpenCL::Accelerator>();

                       /// Destructor.
                       virtual ~MemoryBackend<tags::OpenCL::Accelerator>();

                       virtual void * upload(void * memid, void * address, unsigned long bytes);

                       virtual void download(void * memid, void * address, unsigned long bytes);

                       virtual void * alloc(void * memid, void * address, unsigned long bytes);

                       virtual void free(void * memid);

                       virtual void copy(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                               void * dest_address, unsigned long bytes);

                       virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

                       virtual bool knows(void * memid, void * address);
               };*/
}
#endif
