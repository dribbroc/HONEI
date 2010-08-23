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
#ifndef TOPOLOGY_GUARD_HH
#define TOPOLOGY_GUARD_HH 1

#include <honei/util/attributes.hh>
#include <honei/util/instantiation_policy.hh>

namespace honei
{
    namespace mc
    {
        class Topology :
             public InstantiationPolicy<Topology, Singleton>
        {
            private:
                /// \name Private members
                /// \{

                /// Number of LOGICAL PUs (hardware threads)
                unsigned _num_lpus;

                /// Number of NUMA nodes
                unsigned _num_nodes;

                /// Number of PUs per node
                unsigned _lpus_per_node;

                /// Array that maps sched_ids to NUMA nodes
                unsigned * cpu_to_node;

                /// Lowest scheduler ID of NUMA nodes
                unsigned * range_min;

                /// Highest scheduler ID of NUMA nodes
                unsigned * range_max;

#if defined(__i386__) || defined(__x86_64__)
                /// Processor vendor as specified by the corresponding enumeration
                unsigned _vendor;

                /// Return the number of PUs per physical processor package
                unsigned _num_cores;

                /// Return the number of physical processor packages (num_lpus / num_cores)
                unsigned _num_cpus;

                /// Return the number of hardware threads per processor core (usually 1 or 2)
                unsigned _ht_factor;

#endif
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

#if defined(__i386__) || defined(__x86_64__)

                /// Return the number of PUs per physical processor package
                unsigned num_cores() const;

                /// Return the number of physical processor packages (num_lpus / num_cores)
                unsigned num_cpus() const;

                /// Return the number of hardware threads per processor core (usually 1 or 2)
                unsigned ht_factor() const;
#endif
                /// \}
        };
    }
}
#endif
