/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/backends/memory/memory_backend.hh>
#include <honei/util/instantiation_policy-impl.hh>

namespace honei
{
    template class InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>;
    template class MemoryBackend<tags::CPU>;

    template class InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>;
    template class MemoryBackend<tags::GPU::CUDA>;
}
