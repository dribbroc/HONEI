/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 - 2012 Sven Mallach <mallach@honei.org>
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

#pragma once
#ifndef MULTICORE_GUARD_X86_SPEC_HH
#define MULTICORE_GUARD_X86_SPEC_HH 1

#include <honei/backends/multicore/topology.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>

#include <stdint.h>
#include <unistd.h>

#include <iostream>

namespace honei
{
    namespace mc
    {
        struct x86_unit
        {
            int sched_id;
            int apic_id;
            int topo_id[3]; // 0 = SMT, 1 = CORE, 2 = PACKAGE
        };

#if defined(__i386__)
        inline void cpuid(int CPUInfo[4], int infoType, int ecx_init = 0)
        {
            uint32_t eax, ebx, ecx, edx;
            __asm__ (
                "push %%ebx\n\t"
                "cpuid\n\t"
                "movl %%ebx, %0\n\t"
                "pop %%ebx\n\t"
                : "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
                : "a" (infoType), "c" (ecx_init)
                : "cc"
            );

            CPUInfo[0] = eax;
            CPUInfo[1] = ebx;
            CPUInfo[2] = ecx;
            CPUInfo[3] = edx;
        }

#elif defined(__x86_64__)
        inline void cpuid(int CPUInfo[4], int infoType, int rcx_init = 0)
        {
            uint64_t _rax, _rbx, _rcx, _rdx;
            __asm__ (
                "cpuid\n\t"
                : "=a" (_rax), "=b" (_rbx), "=c" (_rcx), "=d" (_rdx)
                : "a" (infoType), "c" (rcx_init)
                : "cc"
            );

            CPUInfo[0] = static_cast<int>(_rax);
            CPUInfo[1] = static_cast<int>(_rbx);
            CPUInfo[2] = static_cast<int>(_rcx);
            CPUInfo[3] = static_cast<int>(_rdx);
        }
#endif

        inline int next_power_of_two(int src)
        {
            int power = 1;

            while (power < src)
                power <<= 1;

            return power;
        }

        inline int log2(int src)
        {
            int power = 1, count = 0;

            while (power < src)
            {
                ++count;
                power <<= 1;
            }

            return count;
        }

#if defined(__i386__) || defined(__x86_64__)

        // Retrieve the unique APIC processor id of the currently executing core
        template <int> int retrieve_apic_id();

        template <> int retrieve_apic_id<x86_intel>()
        {
            int CPUInfo[4];
            cpuid(CPUInfo, 0x0);

            // check if extended topology leaf (0BH) is available
            if (CPUInfo[0] >= 0xB)
            {
                cpuid(CPUInfo, 0xB);

                if ((CPUInfo[0] & 0x1F) != 1 || CPUInfo[0] == 0 || CPUInfo[1] == 0 || (CPUInfo[2] & 0xFF) != 0)
                    std::cout << "ERROR!" << std::endl;

                int ecx = CPUInfo[2];
                ecx &= 0xFFFF;
                ecx >>= 8;

                if (ecx == 0)
                    std::cout << "ERROR! INVALID LEVEL!" << std::endl;

                std::cout << "EDX = " << CPUInfo[3] << std::endl;

                int result = CPUInfo[3];

                return result; // x2APIC id is in EDX
            }
            else if (CPUInfo[0] >= 0x1)
            {
                cpuid(CPUInfo, 0x1);

                int result = (CPUInfo[1] >> 24); // Read EBX[31:24]
                return result;
            }
            else
            {
                return -1;
            }
        }

        template <> int retrieve_apic_id<x86_amd>()
        {
            return 0;
        }

        // Retrieve the topology for a given unit using leaf OBH
        template <int> x86_unit * retrieve_topology(unsigned sched);

        /*
        template <> x86_unit * retrieve_topology<x86_amd>(unsigned sched)
        {
            int CPUInfo1[4];
            cpuid(CPUInfo1, 0x1);

            int CPUInfo88[4];
            cpuid(CPUInfo88, 0x80000008);

            int lapic = CPUInfo1[1] >> 24;
            int lpc = (CPUInfo1[1] >> 16) & 0xFF;

            int htt = (CPUInfo1[3] >> 28) & 0x1;
            int cmp_legacy = (CPUInfo1[2] >> 1) & 0x1;

            // ECX[15:12] !! Error in CPUID spec document
            int apicIdCoreSize = (CPUInfo88[3] >> 12) & 0x1F;

            if (apicIdCoreSize > 0) // Extended Method available
            {
                int mnc = 1 << apicIdCoreSize;
                int nc = 1 + (CPUInfo88[2] & 0xFF);


            }
            else // legacy method
            {
                int nc = 1 + (CPUInfo88[2] & 0xFF);
                int mnc = nc;
            }
        }
        */

        template <> x86_unit * retrieve_topology<x86_intel>(unsigned sched)
        {
            x86_unit * unit = new x86_unit;

            unit->sched_id = sched;

            int CPUInfo[4];
            cpuid(CPUInfo, 0x0);

            // check if extended topology leaf (0BH) is available
            if (CPUInfo[0] >= 0xB)
            {
                int smt_mask_width(-1), smt_select_mask, core_mask_width(-1);

                cpuid(CPUInfo, 0xB, 0);

                int ecx = CPUInfo[2];
                int ecx_low16 = ecx & 0xFFFF;
                int ecx_15_8 = ecx_low16 >> 8;

                int edx = CPUInfo[3] & 0xFFFFFFFF;
                unit->apic_id = edx;

                if (ecx_15_8 == 1) // ECX[15:8] must be 1
                {
                    smt_mask_width = CPUInfo[0] & 0x1F; // EAX[4:0]

                    if (smt_mask_width != 1)
                        std::cout << "ERROR WHILE RETRIEVING SMT_ID" << std::endl;

                    smt_select_mask = ~((-1) << smt_mask_width);
                    std::cout << "select mask: " << smt_select_mask << std::endl;


                    unit->topo_id[0] = unit->apic_id & smt_select_mask;
                    std::cout << "apic id is : " << unit->apic_id << "  -  got SMT id " << unit->topo_id[0] << std::endl;
                }
                else
                {
                    std::cout << "ERROR WHILE RETRIEVING SMT_ID" << std::endl;
                }

                cpuid(CPUInfo, 0xB, 1);

                ecx = CPUInfo[2];
                ecx_low16 = ecx & 0xFFFF;
                ecx_15_8 = ecx_low16 >> 8;

                if (ecx_15_8 == 2 && smt_mask_width >= 0) // ECX[15:8] must be 2
                {
                    core_mask_width = CPUInfo[0] & 0x1F; // EAX[4:0]
                    std::cout << "core mask width: " << core_mask_width << std::endl;
                    int core_only_select_mask = (~((-1) << core_mask_width)) ^ smt_select_mask;
                    unit->topo_id[1] = (unit->apic_id & core_only_select_mask) >> smt_mask_width;
                    std::cout << "apic id is : " << unit->apic_id << "  -  got core id " << unit->topo_id[1] << std::endl;
                }

                int package_select_mask = ((-1) << core_mask_width);
                unit->topo_id[2] = (unit->apic_id & package_select_mask) >> core_mask_width;
            }
            else if (CPUInfo[0] >= 0x1)
            {
                // Retrieve the topology for a given unit using legacy method
                // Attention: Not working correctly for all processors though
                // implemented in compliance with INTEL manuals.

                unit->apic_id = retrieve_apic_id<x86_intel>();

                cpuid(CPUInfo, 0x1);
                int max_apics = (CPUInfo[1] >> 16) & 0xFF; // EBX[23:16]

                cpuid(CPUInfo, 0x4);
                int max_cores = (CPUInfo[0] >> 26); // EAX[31:26]

                int smt_mask_width = log2(next_power_of_two(max_apics) / (1 + max_cores));
                int smt_select_mask = ~((-1) << smt_mask_width);
                unit->topo_id[0] = unit->apic_id & smt_select_mask;

                int core_mask_width = log2(1 + max_cores);
                int core_select_mask = (~((-1) << (core_mask_width + smt_mask_width))) ^ smt_select_mask;
                unit->topo_id[1] = (unit->apic_id & core_select_mask) >> smt_mask_width;

                int package_mask_width = core_mask_width + smt_mask_width;
                int package_select_mask = ((-1) << package_mask_width);
                unit->topo_id[2] = (unit->apic_id & package_select_mask) >> package_mask_width;
            }

            return unit;
        }


#endif
    }
}
#endif
