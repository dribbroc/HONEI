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

#include <honei/backends/multicore/topology.hh>
#include <honei/util/private_implementation_pattern-impl.hh>

#include <stdint.h>
#include <sys/syscall.h>

namespace honei
{
    namespace mc
    {
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
    }

    using namespace mc;

    template <> struct Implementation<mc::Topology>
    {
        /// Number of LOGICAL PUs (hardware threads)
        unsigned num_lpus;

        /// Number of LOGICAL PUs per physical processor package
        unsigned num_cores;

        /// Number of physical processor packages (num_lpus / num_cores)
        unsigned num_cpus;

        /// Number of hardware threads per processor core (usually 1 or 2)
        unsigned ht_factor;

#if defined(__i386__) || defined(__x86_64__)

        /// Enumeration to distinguish the main x86 processor vendors
        enum x86_Vendor
        {
            UNDEFINED,
            INTEL,
            AMD
        };

        /// Processor vendor as specified by the corresponding enumeration
        int vendor;

        /// Read processor specific information using cpuid instruction
        void init_x86()
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
                    num_cores = 1 + (CPUInfo[0] >> 26); // doubled if HTT is enabled!
                    cpuid(CPUInfo, 0x1);

                    if (((CPUInfo[3] >> 28) & 0x1) == 1)
                        ht_factor = 2;
                    else
                        ht_factor = 1;
                }
            }
            else if (vendor == AMD)
            {
                cpuid(CPUInfo, 0x80000000);
                if (CPUInfo[0] > int(0x80000008))
                {
                    cpuid(CPUInfo, 0x80000008);
                    {
                        num_cores = 1 + (CPUInfo[2] & 0xFF);
                    }

                    ht_factor = 1; // ToDo: Remove hardcoded numbers
                }
            }
        }
#endif

        Implementation()
        {
#if defined(__i386__) || defined(__x86_64__)
            init_x86();
#endif

#if defined linux
            num_lpus = sysconf(_SC_NPROCESSORS_CONF);

            if (vendor == UNDEFINED)
                num_cores = sysconf(_SC_NPROCESSORS_CONF);

#else
            num_lpus = 1;
            num_cores = 1; // ToDo: Remove hardcoded numbers
#endif
            num_cpus = num_lpus / num_cores;
        }
    };
}

using namespace honei::mc;

Topology::Topology() :
    PrivateImplementationPattern<Topology, Single>(new Implementation<Topology>())
{
}

unsigned Topology::num_lpus()
{
    return _imp->num_lpus;
}

#if defined(__i386__) || defined(__x86_64__)

unsigned Topology::num_cores()
{
    return _imp->num_cores;
}

unsigned Topology::num_cpus()
{
    return _imp->num_cpus;
}

unsigned Topology::ht_factor()
{
    return _imp->ht_factor;
}

#endif
