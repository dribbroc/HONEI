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
#include <honei/util/mutex.hh>
#include <honei/util/condition_variable.hh>
#include <honei/backends/memory/memory_backend.hh>

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
            Mutex * const mutex;
            ConditionVariable * access_finished;

            std::map<unsigned long, MemoryBlock> blocks;

        public:
            friend class InstantiationPolicy<MemoryArbiter, Singleton>;

            MemoryArbiter() :
                mutex(new Mutex),
                access_finished(new ConditionVariable)
            {
            }

            ~MemoryArbiter()
            {
                delete access_finished;
                delete mutex;
            }

        template<typename Tag_>
            void * read(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving read lock:");
                {
                    Lock l(*mutex);
                    std::map<unsigned long, MemoryBlock>::iterator i(blocks.find(memid));
                    while (i != blocks.end() && i->second.write_count != 0)
                    {
                        access_finished->wait(*mutex);
                        i = blocks.find(memid);
                    }
                    if (i == blocks.end())
                    {
                        MemoryBlock new_block(tags::tv_none);
                        new_block.read_count++;
                        new_block.readers.insert(tags::TagValue(Tag_::memory_value));
                        blocks.insert(std::pair<unsigned long, MemoryBlock>(memid, new_block));
                    }
                    else
                    {
                        if (i->second.writer != tags::tv_none && i->second.writer != Tag_::memory_value)
                        {
                            switch(i->second.writer)
                            {
                                case tags::tv_cpu:
                                    MemoryBackend<tags::CPU>::instance()->download(memid, address, bytes);
                                    break;
#ifdef HONEI_CUDA
                                case tags::tv_gpu_cuda:
                                    MemoryBackend<tags::GPU::CUDA>::instance()->download(memid, address, bytes);
                                    break;
#endif
                                default:
                                    throw InternalError("MemoryArbiter::read tag not found!");
                            }
                        }
                        i->second.read_count++;
                        i->second.readers.insert(tags::TagValue(Tag_::memory_value));
                    }
                }
                return MemoryBackend<Tag_>::instance()->upload(memid, address, bytes);
            }

        template<typename Tag_>
            void * write(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When retrieving write lock:");
                {
                    Lock l(*mutex);
                    std::map<unsigned long, MemoryBlock>::iterator i(blocks.find(memid));
                    while (i != blocks.end() && (i->second.write_count != 0 || i->second.read_count != 0))
                    {
                        access_finished->wait(*mutex);
                        i = blocks.find(memid);
                    }
                    if (i == blocks.end())
                    {
                        MemoryBlock new_block(Tag_::memory_value);
                        new_block.write_count++;
                        new_block.readers.insert(tags::TagValue(Tag_::memory_value));
                        blocks.insert(std::pair<unsigned long, MemoryBlock>(memid, new_block));
                    }
                    else
                    {
                        if (i->second.writer != tags::tv_none && i->second.writer != Tag_::memory_value)
                        {
                            switch(i->second.writer)
                            {
                                case tags::tv_cpu:
                                    MemoryBackend<tags::CPU>::instance()->download(memid, address, bytes);
                                    break;
#ifdef HONEI_CUDA
                                case tags::tv_gpu_cuda:
                                    MemoryBackend<tags::GPU::CUDA>::instance()->download(memid, address, bytes);
                                    break;
#endif
                                default:
                                    throw InternalError("MemoryArbiter::write writer not found!");
                            }
                        }
                        for (std::set<tags::TagValue>::iterator j(i->second.readers.begin()), j_end(i->second.readers.end()) ;
                                j != j_end ; ++j)
                        {
                            if (*j != Tag_::memory_value)
                            {
                                switch(*j)
                                {
                                    case tags::tv_cpu:
                                        MemoryBackend<tags::CPU>::instance()->reset(memid, address, bytes);
                                        break;
#ifdef HONEI_CUDA
                                    case tags::tv_gpu_cuda:
                                        MemoryBackend<tags::GPU::CUDA>::instance()->reset(memid, address, bytes);
                                        break;
#endif
                                    default:
                                        throw InternalError("MemoryArbiter::write reader not found!");
                                }
                            }

                        }
                        i->second.readers.clear();
                        i->second.writer= Tag_::memory_value;
                        i->second.write_count++;
                        i->second.readers.insert(tags::TagValue(Tag_::memory_value));
                    }
                }
                return MemoryBackend<Tag_>::instance()->upload(memid, address, bytes);
            }

        template<typename Tag_>
            void release_read(unsigned long memid)
            {
                CONTEXT("When releasing read lock:");
                Lock l(*mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(blocks.find(memid));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::release_read MemoryBlock not found!");
                }
                else
                {
                    if (i->second.read_count == 0)
                    {
                        throw InternalError("MemoryArbiter::release_read: read counter is already zero!");
                    }
                    /*else if((temp.writer == tags::tv_cpu || temp.writer == tags::tv_none) && temp.read_count == 1 && temp.write_count == 0)
                    {
                        // delete block
                    }*/
                    else
                    {
                        i->second.read_count--;
                    }
                    access_finished->broadcast();
                }
            }

        template<typename Tag_>
            void release_write(unsigned long memid)
            {
                CONTEXT("When releasing write lock:");
                Lock l(*mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(blocks.find(memid));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::release_write MemoryBlock not found!");
                }
                else
                {
                    if (i->second.write_count == 0)
                    {
                        throw InternalError("MemoryArbiter::release_write: write counter is already zero!");
                    }
                    /*else if((temp.writer == tags::tv_cpu || temp.writer == tags::tv_none) && temp.read_count == 0 && temp.write_count == 1)
                    {
                        // delete block
                    }*/
                    else
                    {
                        i->second.write_count--;
                    }
                    access_finished->broadcast();
                }
            }

        template <typename Tag_>
            void try_download(unsigned long memid, void * address, unsigned long bytes)
            {
                CONTEXT("When trying to download unneeded data:");
                Lock l(*mutex);
                std::map<unsigned long, MemoryBlock>::iterator i(blocks.find(memid));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::try_download MemoryBlock not found!");
                }
                else if (i->second.write_count == 0 && i->second.read_count == 0)
                {
                    MemoryBackend<Tag_>::instance->download(memid, address, bytes);
                    blocks.erase(i);
                }
            }
    };
}
#endif
