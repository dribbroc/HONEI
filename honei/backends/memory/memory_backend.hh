/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MEMORY_GUARD_MEMORY_BACKEND_HH
#define MEMORY_GUARD_MEMORY_BACKEND_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/tags.hh>

namespace honei
{

    template<typename Tag_>
    class MemoryBackend
        //public InstantiationPolicy<MemoryBackend<Tag_>, Singleton>
    {
    };

    template<>
    class MemoryBackend<tags::CPU> :
        public InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>
    {
    public:
        unsigned long upload(unsigned long host_address, unsigned long bytes)
        {
            return host_address;
        }

        void download(unsigned long host_address, unsigned long bytes)
        {
        }

    };

    template<>
    class MemoryBackend<tags::GPU::CUDA> :
        public InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>
    {
    public:
        unsigned long upload(unsigned long host_address, unsigned long bytes);

        void download(unsigned long host_address, unsigned long bytes);

    };
}
#endif
