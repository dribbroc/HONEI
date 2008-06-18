/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

namespace honei
{
    struct MemoryBlock
    {
        public:

            MemoryBlock(unsigned long m, tags::TagValue w) :
                memid(m),
                writer(w)
            {
                read_count = 0;
                write_count = 0;
            }

            unsigned long memid;

            unsigned read_count;
            unsigned write_count;
            tags::TagValue writer;
            std::set<tags::TagValue> readers;
    };

    struct MBCompare
    {
            bool operator() (const MemoryBlock &  left, const MemoryBlock &  right) const
            {
                    return left.memid < right.memid;
            }
    };

    class MemoryArbiter :
        public InstantiationPolicy<MemoryArbiter, Singleton>
    {
        private:
            Mutex * const mutex;
            ConditionVariable * access_finished;

            std::set<MemoryBlock, MBCompare> blocks;

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
            unsigned long read(unsigned long memid, unsigned long address, unsigned long bytes)
            {
                CONTEXT("When retrieving read lock:");
                {
                    Lock l(*mutex);
                    MemoryBlock new_block(memid, tags::tv_none);
                    std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                    while (i != blocks.end() && i->write_count != 0)
                    {
                        access_finished->wait(*mutex);
                        i = blocks.find(new_block);
                    }
                    if (i == blocks.end())
                    {
                        new_block.read_count++;
                        new_block.readers.insert(tags::TagValue(Tag_::tag_value));
                        blocks.insert(new_block);
                    }
                    else
                    {
                        MemoryBlock temp(*i);
                        blocks.erase(i);
                        if (temp.writer != tags::tv_none && temp.writer != Tag_::tag_value)
                        {
                            switch(temp.writer)
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
                        temp.read_count++;
                        new_block.readers.insert(tags::TagValue(Tag_::tag_value));
                        blocks.insert(temp);
                    }
                }
                return MemoryBackend<Tag_>::instance()->upload(memid, address, bytes);
            }

        template<typename Tag_>
            unsigned long write(unsigned long memid, unsigned long address, unsigned long bytes)
            {
                CONTEXT("When retrieving write lock:");
                {
                    Lock l(*mutex);
                    MemoryBlock new_block(memid, Tag_::tag_value);
                    std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                    while (i != blocks.end() && (i->write_count != 0 || i->read_count != 0))
                    {
                        access_finished->wait(*mutex);
                        i = blocks.find(new_block);
                    }
                    if (i == blocks.end())
                    {
                        new_block.write_count++;
                        new_block.readers.insert(tags::TagValue(Tag_::tag_value));
                        blocks.insert(new_block);
                    }
                    else
                    {
                        MemoryBlock temp(*i);
                        blocks.erase(i);
                        if (temp.writer != tags::tv_none && temp.writer != Tag_::tag_value)
                        {
                            switch(temp.writer)
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
                        for (std::set<tags::TagValue>::iterator j(temp.readers.begin()), j_end(temp.readers.end()) ;
                                j != j_end ; ++j)
                        {
                            if (*j != Tag_::tag_value)
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
                        temp.writer= Tag_::tag_value;
                        temp.write_count++;
                        temp.readers.insert(tags::TagValue(Tag_::tag_value));
                        blocks.insert(temp);
                    }
                }
                return MemoryBackend<Tag_>::instance()->upload(memid, address, bytes);
            }

        template<typename Tag_>
            void release_read(unsigned long memid)
            {
                CONTEXT("When releasing read lock:");
                Lock l(*mutex);
                MemoryBlock new_block(memid, tags::tv_none);
                std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::release_read MemoryBlock not found!");
                }
                else
                {
                    MemoryBlock temp(*i);
                    blocks.erase(i);
                    if (temp.read_count == 0)
                    {
                        throw InternalError("MemoryArbiter::release_read: read counter is already zero!");
                    }
                    /*else if((temp.writer == tags::tv_cpu || temp.writer == tags::tv_none) && temp.read_count == 1 && temp.write_count == 0)
                    {
                        // leave block deleted
                    }*/
                    else
                    {
                        temp.read_count--;
                        blocks.insert(temp);
                    }
                    access_finished->broadcast();
                }
            }

        template<typename Tag_>
            void release_write(unsigned long memid)
            {
                CONTEXT("When releasing write lock:");
                Lock l(*mutex);
                MemoryBlock new_block(memid, Tag_::tag_value);
                std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::release_write MemoryBlock not found!");
                }
                else
                {
                    MemoryBlock temp(*i);
                    blocks.erase(i);
                    if (temp.write_count == 0)
                    {
                        throw InternalError("MemoryArbiter::release_write: write counter is already zero!");
                    }
                    /*else if((temp.writer == tags::tv_cpu || temp.writer == tags::tv_none) && temp.read_count == 0 && temp.write_count == 1)
                    {
                        // leave block deleted
                    }*/
                    else
                    {
                        temp.write_count--;
                        blocks.insert(temp);
                    }
                    access_finished->broadcast();
                }
            }

        template <typename Tag_>
            void try_download(unsigned long memid, unsigned long address, unsigned long bytes)
            {
                CONTEXT("When trying to download unneeded data:");
                Lock l(*mutex);
                MemoryBlock new_block(memid, Tag_::tag_value);
                std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                if(i == blocks.end())
                {
                    throw InternalError("MemoryArbiter::try_download MemoryBlock not found!");
                }
                else if (i->write_count == 0 && i->read_count == 0)
                {
                    MemoryBackend<Tag_>::instance->download(memid, address, bytes);
                    blocks.erase(i);
                }
            }
    };
}
#endif
