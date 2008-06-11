/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/backends/memory/memory_arbiter.hh>
#include <honei/util/instantiation_policy-impl.hh>

namespace honei
{
    template class InstantiationPolicy<MemoryArbiter, Singleton>;
    class MemoryArbiter;
}
