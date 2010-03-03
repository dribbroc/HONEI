/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/numainfo.hh>
#include <honei/backends/multicore/topology.hh>
#include <honei/util/instantiation_policy-impl.hh>

#include <limits>
#include <sys/syscall.h>
#include <unistd.h>

using namespace honei;
using namespace honei::mc;

template class InstantiationPolicy<Topology, Singleton>;

// ToDo: May be put in an own header for arch-specific code
#if defined(__i386__) || defined(__x86_64__)
void cpuid(int CPUInfo[4], int infoType, int ecx_init = 0)
{
    uint32_t eax, ebx, ecx, edx;
    __asm__ (
        "pushl %%ebx\n\t"
        "cpuid\n\t"
        "movl %%ebx, %1\n\t"
        "popl %%ebx\n\t"
        : "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (infoType), "c" (ecx_init)
        : "cc"
    );

    CPUInfo[0] = eax;
    CPUInfo[1] = ebx;
    CPUInfo[2] = ecx;
    CPUInfo[3] = edx;
}
#endif

Topology::Topology()
{
#if defined(__i386__) || defined(__x86_64__)
        init_x86();
#endif

#if defined linux
        _num_lpus = sysconf(_SC_NPROCESSORS_CONF);

        _num_nodes = intern::num_nodes();

        if (_num_nodes == 1)
        {
            cpu_to_node = new unsigned[_num_lpus];
            for (unsigned i(0) ; i < _num_lpus ; ++i)
                cpu_to_node[i] = 0;

            range_min = new unsigned[1];
            range_max = new unsigned[1];

            range_min[0] = 0;
            range_max[0] = _num_lpus - 1;
        }
        else
        {
            cpu_to_node = intern::cpu_to_node_array(_num_nodes, _num_lpus);

            for (unsigned i(0) ; i < _num_nodes ; ++i)
            {
                range_min[i] = std::numeric_limits<unsigned>::max();
                range_max[i] = std::numeric_limits<unsigned>::min();
            }

            for (unsigned i(0) ; i < _num_lpus ; ++i)
            {
                if (i < range_min[cpu_to_node[i]])
                    range_min[cpu_to_node[i]] = i;

                if (i > range_min[cpu_to_node[i]])
                    range_max[cpu_to_node[i]] = i;
            }
        }

        if (vendor == UNDEFINED)
            _num_cores = sysconf(_SC_NPROCESSORS_CONF);

#else
        _num_lpus = 1;
        _num_cores = 1; // ToDo: Remove hardcoded numbers
#endif
        _num_cpus = _num_lpus / _num_cores;

}
Topology::~Topology()
{
    delete[] cpu_to_node;
    delete[] range_min;
    delete[] range_max;
}

#if defined(__i386__) || defined(__x86_64__)
void Topology::init_x86()
{
    int CPUInfo[4];

    cpuid(CPUInfo, 0x0);
    if (CPUInfo[0] > 0)
    {
        if (CPUInfo[2] == 1818588270)
            vendor = INTEL;
        else if (CPUInfo[2] == 1145913699)
            vendor = AMD;
        else
            vendor = UNDEFINED;
    }

    if (vendor == INTEL)
    {
        cpuid(CPUInfo, 0x4, 0x0);
        if (CPUInfo[0] > 0)
        {
            _num_cores = 1 + (CPUInfo[0] >> 26); // doubled if HTT is enabled!
            cpuid(CPUInfo, 0x1);

            if (((CPUInfo[3] >> 28) & 0x1) == 1)
                _ht_factor = 2;
            else
                _ht_factor = 1;
        }
    }
    else if (vendor == AMD)
    {
        cpuid(CPUInfo, 0x80000000);
        if (CPUInfo[0] > int(0x80000008))
        {
            cpuid(CPUInfo, 0x80000008);
            {
                _num_cores = 1 + (CPUInfo[2] & 0xFF);
            }

            _ht_factor = 1; // ToDo: Remove hardcoded numbers
        }
    }
}
#endif

unsigned Topology::num_lpus() const
{
    return _num_lpus;
}

unsigned Topology::num_nodes() const
{
    return _num_nodes;
}

unsigned Topology::node_min(unsigned node) const
{
    return range_min[node];
}

unsigned Topology::node_max(unsigned node) const
{
    return range_max[node];
}

unsigned Topology::get_node(unsigned lpu) const
{
    return cpu_to_node[lpu];
}

#if defined(__i386__) || defined(__x86_64__)

unsigned Topology::num_cores() const
{
    return _num_cores;
}

unsigned Topology::num_cpus() const
{
    return _num_cpus;
}

unsigned Topology::ht_factor() const
{
    return _ht_factor;
}

#endif
