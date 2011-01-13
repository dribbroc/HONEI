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

#pragma once
#ifndef UTIL_GUARD_MEMORY_ARBITER_IMPL_HH
#define UTIL_GUARD_MEMORY_ARBITER_IMPL_HH 1

#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/tags.hh>
#include <honei/util/exception.hh>
#include <honei/util/assertion.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>
#include <honei/util/mutex.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/memory_backend_base.hh>
#include <honei/util/type_traits.hh>

#include <cstring>
#include <set>
#include <map>
#include <string>

namespace honei
{
    template <> struct Implementation<MemoryArbiter>
    {
        /// Our mutex.
        Mutex * const _mutex;

        /// Our one read/write lock has faded condition variable.
        ConditionVariable * _access_finished;

        /// Our map of all memory backends
        std::map<tags::TagValue, MemoryBackendBase *> _backends;

        /// Logical representation of a used chunk of memory.
        struct MemoryBlock
        {
            public:
                MemoryBlock(tags::TagValue w) :
                    writer(w)
            {
                read_count = 0;
                write_count = 0;
            }

                unsigned read_count;
                unsigned write_count;
                tags::TagValue writer;
                std::set<tags::TagValue> readers;
        };

        /// Our map of all memory blocks.
        std::map<void *, MemoryBlock> _blocks;

        friend class MemoryBackendRegistrator;

        /// Constructor
        Implementation() :
            _mutex(new Mutex),
            _access_finished(new ConditionVariable)
        {
        }

        /// Destructor
        ~Implementation()
        {
            delete _access_finished;
            delete _mutex;
        }

        void register_address(void * memid)
        {
            CONTEXT("When registering memory block...");
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            if (i != _blocks.end())
            {
                throw InternalError("MemoryArbiter: Duplicate Memory Block!");
            }
            else
            {
                MemoryBlock new_block(tags::tv_none);
                _blocks.insert(std::pair<void *, MemoryBlock>(memid, new_block));
            }
        }

        void remove_address(void * memid)
        {
            CONTEXT("When removing memory block...");
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            if (i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                ASSERT(i->second.read_count == 0 , "Deleting MemoryBlock " + stringify(memid) + " that is still under " + stringify(i->second.read_count) + " read access!");
                ASSERT(i->second.write_count == 0, "Deleting MemoryBlock " + stringify(memid) + " that is still under " + stringify(i->second.write_count) + " write access!");
                // Delete the deprecated memory block in all relevant memory backends
                for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                        j != j_end ; ++j)
                {
                    _backends[*j]->free(memid);
                }
                _blocks.erase(i);
            }
        }

        void copy(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
                void * dest_address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator src_i(_blocks.find(src_id));
            // Wait until no one writes on our src memory block
            while (src_i != _blocks.end() && (src_i->second.write_count != 0))
            {
                _access_finished->wait(*_mutex);
                src_i = _blocks.find(src_id);
            }
            std::map<void *, MemoryBlock>::iterator dest_i(_blocks.find(dest_id));
            // Wait until no one reads or writes on our dest memory block
            while (dest_i != _blocks.end() && (dest_i->second.write_count != 0 || dest_i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                dest_i = _blocks.find(dest_id);
            }

            if (src_i == _blocks.end() || dest_i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // If we can copy just in our local device memory
                if ((src_i->second.writer == tags::tv_none || src_i->second.writer == memory) && _backends[memory]->knows(src_id, src_address))
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        if (*j != memory)
                        {
                            _backends[*j]->free(dest_id);
                        }

                    }
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->alloc(dest_id, dest_address, bytes);
                    _backends[memory]->copy(src_id, src_address, dest_id, dest_address, bytes);
                }

                // If our src memory block was changed on any other remote side, write it back to the main memory
                // and upload it to the dest memory.
                else
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        _backends[*j]->free(dest_id);

                    }
                    if (src_i->second.writer != tags::tv_none)
                        _backends[src_i->second.writer]->download(src_id, src_address, bytes);
                    src_i->second.writer = tags::tv_none;
                    std::memcpy((char *)dest_address, (char *)src_address, bytes);
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->upload(dest_id, dest_address, bytes);
                }
            }
        }

        void convert_float_double(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
                void * dest_address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator src_i(_blocks.find(src_id));
            // Wait until no one writes on our src memory block
            while (src_i != _blocks.end() && (src_i->second.write_count != 0))
            {
                _access_finished->wait(*_mutex);
                src_i = _blocks.find(src_id);
            }
            std::map<void *, MemoryBlock>::iterator dest_i(_blocks.find(dest_id));
            // Wait until no one reads or writes on our dest memory block
            while (dest_i != _blocks.end() && (dest_i->second.write_count != 0 || dest_i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                dest_i = _blocks.find(dest_id);
            }

            if (src_i == _blocks.end() || dest_i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // If we can convert just in our local device memory
                if ((src_i->second.writer == tags::tv_none || src_i->second.writer == memory) && _backends[memory]->knows(src_id, src_address))
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        if (*j != memory)
                        {
                            _backends[*j]->free(dest_id);
                        }

                    }
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->alloc(dest_id, dest_address, bytes * 2);
                    _backends[memory]->convert_float_double(src_id, src_address, dest_id, dest_address, bytes);
                }

                // If our src memory block was changed on any other remote side, write it back to the main memory
                // and upload it to the dest memory.
                else
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        _backends[*j]->free(dest_id);

                    }
                    if (src_i->second.writer != tags::tv_none)
                        _backends[src_i->second.writer]->download(src_id, src_address, bytes);
                    src_i->second.writer = tags::tv_none;
                    TypeTraits<float>::convert((double*) dest_address, (float*)src_address, bytes/sizeof(float));
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->upload(dest_id, dest_address, bytes * 2);
                }
            }
        }

        void convert_double_float(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
                void * dest_address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator src_i(_blocks.find(src_id));
            // Wait until no one writes on our src memory block
            while (src_i != _blocks.end() && (src_i->second.write_count != 0))
            {
                _access_finished->wait(*_mutex);
                src_i = _blocks.find(src_id);
            }
            std::map<void *, MemoryBlock>::iterator dest_i(_blocks.find(dest_id));
            // Wait until no one reads or writes on our dest memory block
            while (dest_i != _blocks.end() && (dest_i->second.write_count != 0 || dest_i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                dest_i = _blocks.find(dest_id);
            }

            if (src_i == _blocks.end() || dest_i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // If we can convert just in our local device memory
                if ((src_i->second.writer == tags::tv_none || src_i->second.writer == memory) && _backends[memory]->knows(src_id, src_address))
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        if (*j != memory)
                        {
                            _backends[*j]->free(dest_id);
                        }

                    }
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->alloc(dest_id, dest_address, bytes/2);
                    _backends[memory]->convert_double_float(src_id, src_address, dest_id, dest_address, bytes);
                }

                // If our src memory block was changed on any other remote side, write it back to the main memory
                // and upload it to the dest memory.
                else
                {
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(dest_i->second.readers.begin()), j_end(dest_i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        _backends[*j]->free(dest_id);

                    }
                    if (src_i->second.writer != tags::tv_none)
                        _backends[src_i->second.writer]->download(src_id, src_address, bytes);
                    src_i->second.writer = tags::tv_none;
                    TypeTraits<double>::convert((float*)dest_address, (double*)src_address, bytes/sizeof(double));
                    dest_i->second.readers.clear();
                    dest_i->second.writer = memory;
                    dest_i->second.readers.insert(memory);
                    _backends[memory]->upload(dest_id, dest_address, bytes / 2);
                }
            }
        }

        void fill(tags::TagValue memory, void * memid, void * address, unsigned long bytes, float proto)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            // Wait until no one reads or writes on our memory block
            while (i != _blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                i = _blocks.find(memid);
            }
            if (i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // Delete the deprecated memory block in all relevant memory backends
                for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                        j != j_end ; ++j)
                {
                    if (*j != memory)
                    {
                        _backends[*j]->free(memid);
                    }

                }
                i->second.readers.clear();
                i->second.writer = memory;
                i->second.readers.insert(memory);
                _backends[memory]->alloc(memid, address, bytes);
                _backends[memory]->fill(memid, address, bytes, proto);
            }
        }

        void * read_only(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            // Wait until no one writes on our memory block
            while (i != _blocks.end() && i->second.write_count != 0)
            {
                _access_finished->wait(*_mutex);
                i = _blocks.find(memid);
            }
            if (i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // If our memory block was changed on any other remote side, write it back to the main memory.
                if (i->second.writer != tags::tv_none && i->second.writer != memory)
                {
                    _backends[i->second.writer]->download(memid, address, bytes);
                    i->second.writer = tags::tv_none;
                }
                i->second.read_count++;
                i->second.readers.insert(memory);
            }
            return _backends[memory]->upload(memid, address, bytes);
        }

        void * read_and_write(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            // Wait until no one reads or writes on our memory block
            while (i != _blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                i = _blocks.find(memid);
            }
            if (i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // If our memory block was changed on any other remote side, write it back to the main memory.
                if (i->second.writer != tags::tv_none && i->second.writer != memory)
                {
                    _backends[i->second.writer]->download(memid, address, bytes);
                }
                // Delete the deprecated memory block in all relevant memory backends
                for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                        j != j_end ; ++j)
                {
                    if (*j != memory)
                    {
                        _backends[*j]->free(memid);
                    }

                }
                i->second.readers.clear();
                i->second.writer = memory;
                i->second.write_count++;
                i->second.readers.insert(memory);
            }
            return _backends[memory]->upload(memid, address, bytes);
        }

        void * write_only(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            // Wait until no one reads or writes on our memory block
            while (i != _blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
            {
                _access_finished->wait(*_mutex);
                i = _blocks.find(memid);
            }
            if (i == _blocks.end())
            {
                throw InternalError("MemoryArbiter: Memory Block not found!");
            }
            else
            {
                // Delete the deprecated memory block in all relevant memory backends
                for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                        j != j_end ; ++j)
                {
                    if (*j != memory)
                    {
                        _backends[*j]->free(memid);
                    }

                }
                i->second.readers.clear();
                i->second.writer = memory;
                i->second.write_count++;
                i->second.readers.insert(memory);
            }
            return _backends[memory]->alloc(memid, address, bytes);
        }

        void release_read(void * memid)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            if(i == _blocks.end())
            {
                throw InternalError("MemoryArbiter::release_read MemoryBlock not found!");
            }
            else
            {
                if (i->second.read_count == 0)
                {
                    throw InternalError("MemoryArbiter::release_read: read counter is already zero!");
                }
                else
                {
                    i->second.read_count--;
                }
                _access_finished->broadcast();
            }
        }

        void release_write(void * memid)
        {
            Lock l(*_mutex);
            std::map<void *, MemoryBlock>::iterator i(_blocks.find(memid));
            if(i == _blocks.end())
            {
                throw InternalError("MemoryArbiter::release_write MemoryBlock not found!");
            }
            else
            {
                if (i->second.write_count == 0)
                {
                    throw InternalError("MemoryArbiter::release_write: write counter is already zero!");
                }
                else
                {
                    i->second.write_count--;
                }
                _access_finished->broadcast();
            }
        }
    };

    MemoryArbiter::MemoryArbiter() :
        PrivateImplementationPattern<MemoryArbiter, Shared>(new Implementation<MemoryArbiter>())
    {
    }

    MemoryArbiter::~MemoryArbiter()
    {
    }

    void MemoryArbiter::register_address(void * memid)
    {
        CONTEXT("When registering memory address:");
        _imp->register_address(memid);
    }

    void MemoryArbiter::remove_address(void * memid)
    {
        CONTEXT("When removing memory address:");
        _imp->remove_address(memid);
    }

    void MemoryArbiter::copy(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        _imp->copy(memory, src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryArbiter::convert_float_double(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        _imp->convert_float_double(memory, src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryArbiter::convert_double_float(tags::TagValue memory, void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        _imp->convert_double_float(memory, src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryArbiter::fill(tags::TagValue memory, void * memid, void * address, unsigned long bytes, float proto)
    {
        _imp->fill(memory, memid, address, bytes, proto);
    }

    void * MemoryArbiter::lock(LockMode mode, tags::TagValue memory, void * memid, void * address, unsigned long bytes)
    {
        switch (mode)
        {
            case lm_read_only:
                return read_only(memory, memid, address, bytes);
            case lm_write_only:
                return write_only(memory, memid, address, bytes);
            case lm_read_and_write:
                return read_and_write(memory, memid, address, bytes);
            default:
                throw InternalError("Memory Arbiter: lock mode unknown!");
        }
    }

    void MemoryArbiter::unlock(LockMode mode, void * memid)
    {
        switch (mode)
        {
            case lm_read_only:
                release_read(memid);
                break;
            case lm_write_only:
                release_write(memid);
                break;
            case lm_read_and_write:
                release_write(memid);
                break;
            default:
                throw InternalError("Memory Arbiter: lock mode unknown!");
        }
    }

    void * MemoryArbiter::read_only(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving read_only lock:");
        return _imp->read_only(memory, memid, address, bytes);
    }

    void * MemoryArbiter::read_and_write(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving read_and_write lock:");
        return _imp->read_and_write(memory, memid, address, bytes);
    }

    void * MemoryArbiter::write_only(tags::TagValue memory, void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving write_only lock:");
        return _imp->write_only(memory, memid, address, bytes);
    }

    void MemoryArbiter::release_read(void * memid)
    {
        CONTEXT("When releasing read lock:");
        _imp->release_read(memid);
    }

    void MemoryArbiter::release_write(void * memid)
    {
        CONTEXT("When releasing write lock:");
        _imp->release_write(memid);
    }

    void MemoryArbiter::insert_backend(std::pair<tags::TagValue, MemoryBackendBase *> backend)
    {
        CONTEXT("When adding a new memory backend:");
        Lock l(*this->_imp->_mutex);
        this->_imp->_backends.insert(backend);
    }
}
#endif
