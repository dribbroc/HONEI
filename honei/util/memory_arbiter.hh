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

    class MemoryArbiter :
        public InstantiationPolicy<MemoryArbiter, Singleton>
    {
        private:
            Mutex * const _mutex;
            ConditionVariable * _access_finished;

            std::map<unsigned long, MemoryBlock> _blocks;
            std::map<tags::TagValue, MemoryBackendBase *> _backends;

        public:
            friend class InstantiationPolicy<MemoryArbiter, Singleton>;
            friend class MemoryBackendRegistrator;

            MemoryArbiter() :
                _mutex(new Mutex),
                _access_finished(new ConditionVariable)
            {
            }

            ~MemoryArbiter()
            {
                delete _access_finished;
                delete _mutex;
            }

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
                for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                        j != j_end ; ++j)
                {
                    _backends[*j]->free(memid);
                }
                _blocks.erase(i);
            }
        }


        template<typename Tag_>
            void * read(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving read lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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

        template<typename Tag_>
            void * write(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving write lock:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
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
                    if (i->second.writer != tags::tv_none && i->second.writer != Tag_::memory_value)
                    {
                        _backends[i->second.writer]->download(memid, address, bytes);
                    }
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

        /*template <typename Tag_>
            void try_download(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When trying to download unneeded data:");
                Lock l(*_mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(_blocks.find(memid));
                if(i == _blocks.end())
                {
                    throw InternalError("MemoryArbiter::try_download MemoryBlock not found!");
                }
                else if (i->second.write_count == 0 && i->second.read_count == 0)
                {
                    //MemoryBackend<Tag_>::instance->download(memid, address, bytes);
                    _backends[Tag_::memory_value]->download(memid, address, bytes);
                    _blocks.erase(i);
                }
            }*/
    };

    /**
     * MemoryBackendRegistrator registers a descendant of MemoryBackend with MemoryManager.
     *
     * \ingroup grpmemorymanager
     */
    struct MemoryBackendRegistrator
    {
        /**
         * Constructor.
         *
         * \param v Tag value that the backend is associated with.
         * \param f Singleton-instance pointer to the backend.
         */
        MemoryBackendRegistrator(const tags::TagValue v, MemoryBackendBase * m)
        {
            MemoryArbiter::instance()->_backends.insert(std::make_pair(v, m));
        }
    };
}
#endif
