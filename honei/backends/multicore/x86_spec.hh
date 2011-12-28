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

#pragma once
#ifndef MULTICORE_GUARD_X86_SPEC_HH
#define MULTICORE_GUARD_X86_SPEC_HH 1

#include <honei/util/log.hh>
#include <honei/util/stringify.hh>

#include <stdint.h>
#include <unistd.h>

namespace honei
{
    namespace mc
    {
        /// Enumeration to distinguish the main x86 processor vendors
        enum x86_Vendor
        {
            UNDEFINED,
            INTEL,
            AMD
        };

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

#if defined(__i386__) || defined(__x86_64__)
        /// Read processor specific information using cpuid instruction
        void init_x86(unsigned & vendor, unsigned & num_cores, unsigned & ht_factor)
        {
            CONTEXT("When retrieving x86-specific information from a processor via cpuid:");

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

#ifdef DEBUG
                    std::string msg = "Found INTEL processor(s) with " + stringify(num_cores) + " LPUs and HTT " +
                    (ht_factor == 1 ? "available" : "not available") + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                }
#ifdef DEBUG
                else
                {
                    std::string msg = "Found INTEL processor(s) but could not further evaluate them\n";
                    LOGMESSAGE(lc_backend, msg);
                }
#endif

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
#ifdef DEBUG
                    std::string msg = "Found AMD processor(s) with " + stringify(num_cores) + " LPUs and HTT " +
                    (ht_factor == 1 ? "available" : "not available") + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                }
#ifdef DEBUG
                else
                {
                    std::string msg = "Found AMD processor(s) but could not further evaluate them\n";
                    LOGMESSAGE(lc_backend, msg);
                }
#endif
            }
#ifdef DEBUG
            else
            {
                std::string msg = "Could not determine vendor and further information about available processors\n";
                LOGMESSAGE(lc_backend, msg);
            }
#endif
        }
#endif

#if defined(__i386__) || defined(__x86_64__)
        // Retrieve the unique APIC processor id of the currently executing core
        inline int retrieve_apic_id_11()
        {
            int CPUInfo[4];
            cpuid(CPUInfo, 0x0);

            // check if extended topology leaf (0BH) is available
            if (CPUInfo[0] >= 0xB)
            {
                cpuid(CPUInfo, 0xB);

                if ((CPUInfo[1] & 0xFFFF) != 0)
                    return CPUInfo[3]; // x2APIC id is in EDX
            }

            return -1;
        }

        inline int retrieve_apic_id_1()
        {
            int CPUInfo[4];

            // Read initial APIC id with old legacy method
            cpuid(CPUInfo, 0x0);
            if (CPUInfo[0] >= 0x1)
            {
                cpuid(CPUInfo, 0x1);
                return (CPUInfo[1] >> 24); // Read EBX[31:24]
            }

            return -1;
        }

        // Retrieve the topology for a given unit using leaf OBH
        inline x86_unit * retrieve_topology_11(unsigned sched)
        {
            x86_unit * unit = new x86_unit;

            unit->sched_id = sched;
            unit->apic_id = retrieve_apic_id_11();

            int CPUInfo[4];
            int smt_mask_width, smt_select_mask, core_mask_width;

            cpuid(CPUInfo, 0xB, 0);
            if (((CPUInfo[2] >> 8) & 0xFF) == 1) // ECX[15:8] must be 1
            {
                smt_mask_width = CPUInfo[0] & 0x1F; // EAX[4:0]
                smt_select_mask = ~((-1) << smt_mask_width);
                unit->topo_id[0] = unit->apic_id & smt_select_mask;
            }

            cpuid(CPUInfo, 0xB, 1);
            if (((CPUInfo[2] >> 8) & 0xFF) == 2) // ECX[15:8] must be 1
            {
                core_mask_width = CPUInfo[0] & 0x1F; // EAX[4:0]
                int core_select_mask = ((~(-1)) << core_mask_width) ^ smt_select_mask;
                unit->topo_id[1] = (unit->apic_id & core_select_mask) >> smt_mask_width;
            }

            int package_select_mask = ((-1) << core_mask_width);
            unit->topo_id[2] = (unit->apic_id & package_select_mask) >> core_mask_width;

            return unit;
        }

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

        // Retrieve the topology for a given unit using legacy method
        // Attention: Not working correctly for all processors though
        // implemented in compliance with INTEL manuals.
        inline x86_unit * retrieve_topology_1(unsigned sched)
        {
            x86_unit * unit = new x86_unit;

            unit->sched_id = sched;
            unit->apic_id = retrieve_apic_id_1();

            int CPUInfo[4];

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

            return unit;
        }

#endif
    }
}
#endif
