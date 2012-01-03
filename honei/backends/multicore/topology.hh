/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Sven Mallach <mallach@honei.org>
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
#ifndef TOPOLOGY_GUARD_HH
#define TOPOLOGY_GUARD_HH 1

#include <honei/backends/multicore/lpu.hh>
#include <honei/util/attributes.hh>
#include <honei/util/barrier.hh>
#include <honei/util/instantiation_policy.hh>

namespace honei
{
    namespace mc
    {
        enum Architectures
        {
            unknown = 0,
            x86_intel,
            x86_amd
        };

        template <int> struct TopologyThreadFunction;

        template <> struct TopologyThreadFunction<x86_intel>
        {
            LPU * const lpu;
            int apic_id;
            Barrier * barrier;

            TopologyThreadFunction(LPU * const pu, Barrier * b);

            void operator() ();
        };


        class Topology :
             public InstantiationPolicy<Topology, Singleton>
        {
            private:
                /// \name Private members
                /// \{

                /// Number of LOGICAL PUs (hardware threads)
                unsigned _num_lpus;

                /// Return the number of PUs per physical processor package
                unsigned _num_cores;

                /// Return the number of physical processor packages (num_lpus / num_cores)
                /// Equivalent to the number of nodes detected in case of NUMA
                unsigned _num_cpus;

                /// Processor arch and vendor (if relevant)
                unsigned _arch;

                /// Array with access to LPU data structures
                LPU ** _lpus;

                /// Socket / Node information
                Socket ** _sockets;

#if defined(__i386__) || defined(__x86_64__)

                /// Return whether the processor support simultaneous multithreading
                bool _ht_support;

                /// Return the number of hardware threads per processor core (usually 1 or 2)
                unsigned _ht_factor;

#endif
                /// Determine the architecture of the underlying system
                void determine_arch();

                void enumerate_x86_intel();
                void enumerate_x86_amd();
                void enumerate_numainfo(int num_nodes);

                /// \}

            protected:

                friend class InstantiationPolicy<Topology, Singleton>;

                /// \name Basic Operations
                /// \{

                /// Constructor
                Topology();

                /// \}

            public:
                /// \name Basic Operations
                /// \{

                /// Destructor
                ~Topology();

                /// \}

                /// \name Public members
                /// \{

                LPU ** lpus();
                LPU * lpu(int id);

                Socket ** sockets();

                /// Return the number of logical PUs (hardware-threads)
                unsigned num_lpus() const;

                /// Return the number of NUMA nodes
                unsigned num_nodes() const;

                /// Return the number of PUs per node
                unsigned lpus_per_node() const;

                /// Return the node the lpu belongs to
                unsigned get_node(unsigned lpu) const;

                /// Return the lowest sched_id of a lpu belong to the node
                unsigned node_min(unsigned node) const;

                /// Return the highest sched_id of a lpu belong to the node
                unsigned node_max(unsigned node) const;

                /// Return the node which the main thread is running on
                unsigned main_node() const;

                /// Return the number of PUs per physical processor package
                unsigned num_cores() const;

                /// Return the number of physical processor packages (num_lpus / num_cores)
                unsigned num_cpus() const;

#if defined(__i386__) || defined(__x86_64__)


                /// Return the number of hardware threads per processor core (usually 1 or 2)
                unsigned ht_factor() const;
#endif
                /// \}
        };
    }
}
#endif
