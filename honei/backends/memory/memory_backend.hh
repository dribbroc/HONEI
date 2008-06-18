/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MEMORY_GUARD_MEMORY_BACKEND_HH
#define MEMORY_GUARD_MEMORY_BACKEND_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/tags.hh>

#include <map>

namespace honei
{

    template<typename Tag_>
    class MemoryBackend
    {
    };

    template<>
    class MemoryBackend<tags::CPU> :
        public InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>
    {
    public:
        unsigned long upload(unsigned long memid, unsigned long address, unsigned long size)
        {
            return address;
        }

        void download(unsigned long memid, unsigned long address, unsigned long size)
        {
        }

        void reset(unsigned long memid, unsigned long address, unsigned long size)
        {
        }

    };

    template<>
    class MemoryBackend<tags::GPU::CUDA> :
        public InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>
    {
        private:
            std::map<unsigned long, unsigned long> address_map;

        public:
            unsigned long upload(unsigned long memid, unsigned long address, unsigned long size);

            void download(unsigned long memid, unsigned long address, unsigned long size);

            void reset(unsigned long memid, unsigned long address, unsigned long size);
    };
}
#endif
