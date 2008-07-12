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

#ifndef MEMORY_GUARD_MEMORY_CONTROLLER_IMPL_HH
#define MEMORY_GUARD_MEMORY_CONTROLLER_IMPL_HH 1

#include <honei/util/instantiation_policy.hh>
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

#include <set>
#include <map>

namespace honei
{
    template <> struct Implementation<MemoryArbiter>
    {
        /// Our mutex.
        Mutex * const _mutex;

        /// Our one read/write lock has faded condition variable.
        ConditionVariable * _access_finished;

        /// Our map of all memory blocks.
        std::map<unsigned long, MemoryBlock> _blocks;

        /// Our map of all memory backends
        std::map<tags::TagValue, MemoryBackendBase *> _backends;

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
    };

    MemoryArbiter::MemoryArbiter() :
        PrivateImplementationPattern<MemoryArbiter, Shared>(new Implementation<MemoryArbiter>())
    {
    }

    MemoryArbiter::~MemoryArbiter()
    {
    }

    void MemoryArbiter::add_memblock(unsigned long memid)
    {
        CONTEXT("When adding Memory Block:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        if (i != this->_imp->_blocks.end())
        {
            throw InternalError("MemoryArbiter: Duplicate Memory Block!");
        }
        else
        {
            MemoryBlock new_block(tags::tv_none);
            this->_imp->_blocks.insert(std::pair<unsigned long, MemoryBlock>(memid, new_block));
        }
    }

    void MemoryArbiter::remove_memblock(unsigned long memid)
    {
        CONTEXT("When removing Memory Block:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        if (i == this->_imp->_blocks.end())
        {
            throw InternalError("MemoryArbiter: Memory Block not found!");
        }
        else
        {
            //ASSERT(i->second.read_count == 0 , "Deleting MemoryBlock " + stringify(memid) + " that is still under " + stringify(i->second.read_count) + " read access!");
            //ASSERT(i->second.write_count == 0, "Deleting MemoryBlock " + stringify(memid) + " that is still under " + stringify(i->second.write_count) + " write access!");
            // Delete the deprecated memory block in all relevant memory backends
            for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                    j != j_end ; ++j)
            {
                this->_imp->_backends[*j]->free(memid);
            }
            this->_imp->_blocks.erase(i);
        }
    }

    void * MemoryArbiter::lock(LockMode mode, tags::TagValue memory, unsigned long memid, void * address, unsigned long bytes)
    {
        switch (mode)
        {
            case lm_read_only:
                return read_only(memory, memid, address, bytes);
            case lm_write_only:
                return write_only(memory, memid, address, bytes);
            case lm_read_and_write:
                return read_and_write(memory, memid, address, bytes);
        }
    }

    void MemoryArbiter::unlock(LockMode mode, unsigned long memid)
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
        }
    }

    void * MemoryArbiter::read_only(tags::TagValue memory, unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving read lock:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        // Wait until no one writes on our memory block
        while (i != this->_imp->_blocks.end() && i->second.write_count != 0)
        {
            this->_imp->_access_finished->wait(*this->_imp->_mutex);
            i = this->_imp->_blocks.find(memid);
        }
        if (i == this->_imp->_blocks.end())
        {
            throw InternalError("MemoryArbiter: Memory Block not found!");
        }
        else
        {
            // If our memory block was changed on any other remote side, write it back to the main memory.
            if (i->second.writer != tags::tv_none && i->second.writer != memory)
            {
                this->_imp->_backends[i->second.writer]->download(memid, address, bytes);
                i->second.writer = tags::tv_none;
            }
            i->second.read_count++;
            i->second.readers.insert(memory);
        }
        return this->_imp->_backends[memory]->upload(memid, address, bytes);
    }

    void * MemoryArbiter::read_and_write(tags::TagValue memory, unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving write lock:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        // Wait until no one reads or writes on our memory block
        while (i != this->_imp->_blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
        {
            this->_imp->_access_finished->wait(*this->_imp->_mutex);
            i = this->_imp->_blocks.find(memid);
        }
        if (i == this->_imp->_blocks.end())
        {
            throw InternalError("MemoryArbiter: Memory Block not found!");
        }
        else
        {
            // If our memory block was changed on any other remote side, write it back to the main memory.
            if (i->second.writer != tags::tv_none && i->second.writer != memory)
            {
                this->_imp->_backends[i->second.writer]->download(memid, address, bytes);
            }
            // Delete the deprecated memory block in all relevant memory backends
            for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                    j != j_end ; ++j)
            {
                if (*j != memory)
                {
                    this->_imp->_backends[*j]->free(memid);
                }

            }
            i->second.readers.clear();
            i->second.writer= memory;
            i->second.write_count++;
            i->second.readers.insert(memory);
        }
        return this->_imp->_backends[memory]->upload(memid, address, bytes);
    }

    void * MemoryArbiter::write_only(tags::TagValue memory, unsigned long memid, void * address, unsigned long bytes)
    {
        CONTEXT("When retrieving write-only lock:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        // Wait until no one reads or writes on our memory block
        while (i != this->_imp->_blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
        {
            this->_imp->_access_finished->wait(*this->_imp->_mutex);
            i = this->_imp->_blocks.find(memid);
        }
        if (i == this->_imp->_blocks.end())
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
                    this->_imp->_backends[*j]->free(memid);
                }

            }
            i->second.readers.clear();
            i->second.writer= memory;
            i->second.write_count++;
            i->second.readers.insert(memory);
        }
        return this->_imp->_backends[memory]->alloc(memid, address, bytes);
    }

    void MemoryArbiter::release_read(unsigned long memid)
    {
        CONTEXT("When releasing read lock:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        if(i == this->_imp->_blocks.end())
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
            this->_imp->_access_finished->broadcast();
        }
    }

    void MemoryArbiter::release_write(unsigned long memid)
    {
        CONTEXT("When releasing write lock:");
        Lock l(*this->_imp->_mutex);
        std::map<unsigned long, MemoryBlock>::iterator i(this->_imp->_blocks.find(memid));
        if(i == this->_imp->_blocks.end())
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
            this->_imp->_access_finished->broadcast();
        }
    }

    void MemoryArbiter::insert_backend(std::pair<tags::TagValue, MemoryBackendBase *> backend)
    {
        CONTEXT("When adding a new memory backend:");
        Lock l(*this->_imp->_mutex);
        this->_imp->_backends.insert(backend);
    }


}
#endif
