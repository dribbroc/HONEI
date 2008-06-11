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

            MemoryBlock(unsigned long a, unsigned long s, tags::TagValue v) :
                address(a),
                size(s),
                value(v)
            {
                read_count = 0;
                write_count = 0;
            }

            unsigned long address;
            unsigned long size;

            unsigned read_count;
            unsigned write_count;
            tags::TagValue value;
    };

    struct MBCompare
    {
            bool operator() (const MemoryBlock &  left, const MemoryBlock &  right) const
            {
                //if(left.address != right.address)
                    return left.address < right.address;

                //else
                //    return left.value < right.value;
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


        /// \todo ops blocken bis rdy und downloaden wenn noetig
        /// \todo set ab und zu bereinigen
        template<typename Tag_>
            unsigned long read(unsigned long address, unsigned long bytes)
            {
                Lock l(*mutex);
                MemoryBlock new_block(address, bytes, tags::tv_none);
                std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                while (i != blocks.end() && i->write_count != 0)
                {
                    access_finished->wait(*mutex);
                    i = blocks.find(new_block);
                }
                if (i == blocks.end())
                {
                    new_block.read_count++;
                    blocks.insert(new_block);
                }
                else
                {
                    MemoryBlock temp(*i);
                    blocks.erase(i);
                    if (temp.value != tags::tv_none && temp.value != Tag_::tag_value)
                    {
                        switch(temp.value)
                        {
                            case tags::tv_cpu:
                                MemoryBackend<tags::CPU>::instance()->download(temp.address, temp.size);
                                break;
                            case tags::tv_gpu_cuda:
                        //        MemoryBackend<tags::GPU::CUDA>::instance()->download(temp.address, temp.size);
                                break;
                        }
                    }
                    temp.read_count++;
                    temp.value = tags::tv_none;
                    blocks.insert(temp);
                }
                return MemoryBackend<Tag_>::instance()->upload(address, bytes);
            }

        template<typename Tag_>
            unsigned long write(unsigned long address, unsigned long bytes)
            {
                Lock l(*mutex);
                MemoryBlock new_block(address, bytes, Tag_::tag_value);
                std::set<MemoryBlock, MBCompare>::iterator i(blocks.find(new_block));
                while (i != blocks.end() && (i->write_count != 0 || i->read_count != 0))
                {
                    access_finished->wait(*mutex);
                    i = blocks.find(new_block);
                }
                if (i == blocks.end())
                {
                    new_block.write_count++;
                    blocks.insert(new_block);
                }
                else
                {
                    MemoryBlock temp(*i);
                    blocks.erase(i);
                    if (temp.value != tags::tv_none && temp.value != Tag_::tag_value)
                    {
                        switch(temp.value)
                        {
                            case tags::tv_cpu:
                                MemoryBackend<tags::CPU>::instance()->download(temp.address, temp.size);
                                break;
                            case tags::tv_gpu_cuda:
                          //      MemoryBackend<tags::GPU::CUDA>::instance()->download(temp.address, temp.size);
                                break;
                        }
                    }
                    temp.value= Tag_::tag_value;
                    temp.write_count++;
                    blocks.insert(temp);
                }
                return MemoryBackend<Tag_>::instance()->upload(address, bytes);
            }

        template<typename Tag_>
            void release_read(unsigned long address, unsigned long bytes)
            {
                Lock l(*mutex);
                MemoryBlock new_block(address, bytes, tags::tv_none);
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
                    else
                    {
                    temp.read_count--;
                        blocks.insert(temp);
                    }
                    access_finished->broadcast();
                }
            }

        template<typename Tag_>
            void release_write(unsigned long address, unsigned long bytes)
            {
                Lock l(*mutex);
                MemoryBlock new_block(address, bytes, Tag_::tag_value);
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
                    else
                    {
                        temp.write_count--;
                        blocks.insert(temp);
                    }
                    access_finished->broadcast();
                }
            }
    };
}
#endif
