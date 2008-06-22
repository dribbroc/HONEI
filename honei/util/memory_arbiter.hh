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

#ifndef MEMORY_GUARD_MEMORY_CONTROLLER_HH
#define MEMORY_GUARD_MEMORY_CONTROLLER_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/tags.hh>
#include <honei/util/exception.hh>
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

    /**
     * MemoryArbiter handles read/write locks for all used memory blocks and
     * distributes necessary memory transfers jobs.
     *
     * \ingroup grpmemorymanager
     */
    class MemoryArbiter :
        public InstantiationPolicy<MemoryArbiter, Singleton>
    {
        private:
            /// Our mutex.
            Mutex * const _mutex;

            /// Our one read/write lock has faded condition variable.
            ConditionVariable * _access_finished;

            /// Our map of all memory blocks.
            std::map<unsigned long, MemoryBlock> _blocks;

            /// Our map of all memory backends
            std::map<tags::TagValue, MemoryBackendBase *> _backends;

        public:
            friend class InstantiationPolicy<MemoryArbiter, Singleton>;
            friend class MemoryBackendRegistrator;

            /**
             * Constructor.
             */
            MemoryArbiter() :
                _mutex(new Mutex),
                _access_finished(new ConditionVariable)
            {
            }

            /**
             * Destructor.
             */
            ~MemoryArbiter()
            {
                delete _access_finished;
                delete _mutex;
            }

            /**
             * Add a chunk of memory to the memory block map.
             *
             * \param memid A unique key identifying the added chunk.
             */
            void add_memblock(unsigned long memid)
            {
                CONTEXT("When adding Memory Block:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
                if (i != _blocks.end())
                {
                    throw InternalError("MemoryArbiter: Duplicate Memory Block!");
                }
                else
                {
                    MemoryBlock new_block(tags::tv_none);
                    _blocks.insert(std::pair<unsigned long, MemoryBlock>(memid, new_block));
                }
            }

            /**
             * Remove a chunk of memory from the memory block map.
             *
             * \param memid A unique key identifying the removed chunk.
             */
            void remove_memblock(unsigned long memid)
            {
                CONTEXT("When removing Memory Block:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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
                        _backends[*j]->free(memid);
                    }
                    _blocks.erase(i);
                }
            }

            /**
             * Request a read lock for a specific memory block.
             *
             * \param memid A unique key identifying the requested chunk.
             * \param address The address where our reading will begin.
             * \param bytes The amount of bytes we want to read.
             */
            template<typename Tag_>
            void * read(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving read lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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
                    if (i->second.writer != tags::tv_none && i->second.writer != Tag_::memory_value)
                    {
                        _backends[i->second.writer]->download(memid, address, bytes);
                        i->second.writer = tags::tv_none;
                    }
                    i->second.read_count++;
                    i->second.readers.insert(tags::TagValue(Tag_::memory_value));
                }
                return _backends[Tag_::memory_value]->upload(memid, address, bytes);
            }

            /**
             * Request a write lock for a specific memory block.
             *
             * \param memid A unique key identifying the requested chunk.
             * \param address The address where our writing will begin.
             * \param bytes The amount of bytes we want to write.
             */
            template<typename Tag_>
            void * write(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving write lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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
                    if (i->second.writer != tags::tv_none && i->second.writer != Tag_::memory_value)
                    {
                        _backends[i->second.writer]->download(memid, address, bytes);
                    }
                    // Delete the deprecated memory block in all relevant memory backends
                    for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                            j != j_end ; ++j)
                    {
                        if (*j != Tag_::memory_value)
                        {
                            _backends[*j]->free(memid);
                        }

                    }
                    i->second.readers.clear();
                    i->second.writer= Tag_::memory_value;
                    i->second.write_count++;
                    i->second.readers.insert(tags::TagValue(Tag_::memory_value));
                }
                return _backends[Tag_::memory_value]->upload(memid, address, bytes);
            }

            /**
             * Release a read lock for a specific memory block.
             *
             * \param memid A unique key identifying the released chunk.
             */
            template<typename Tag_>
            void release_read(unsigned long memid)
            {
                CONTEXT("When releasing read lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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

            /**
             * Release a write lock for a specific memory block.
             *
             * \param memid A unique key identifying the released chunk.
             */
            template<typename Tag_>
            void release_write(unsigned long memid)
            {
                CONTEXT("When releasing write lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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

    /**
     * MemoryBackendRegistrator registers a descendant of MemoryBackendBase with MemoryArbiter.
     *
     * \ingroup grpmemorymanager
     */
    struct MemoryBackendRegistrator
    {
        /**
         * Constructor.
         *
         * \param v Memory tag value that the backend is associated with.
         * \param f Singleton-instance pointer to the backend.
         */
        MemoryBackendRegistrator(const tags::TagValue v, MemoryBackendBase * m)
        {
            MemoryArbiter::instance()->_backends.insert(std::make_pair(v, m));
        }
    };
}
#endif
