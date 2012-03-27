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

#include <honei/backends/multicore/numainfo.hh>
#include <honei/backends/multicore/topology.hh>
#if defined(__i386__) || defined(__x86_64__)
#include <honei/backends/multicore/x86_spec.hh>
#endif

#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread.hh>

#include <errno.h>
#include <limits>
#include <sched.h>
#include <sys/syscall.h>

using namespace honei;
using namespace honei::mc;

#if defined(__i386__) || defined(__x86_64__)
TopologyThreadFunction<x86_intel>::TopologyThreadFunction(LPU * const pu, Barrier * b) :
    lpu(pu),
    apic_id(-1),
    barrier(b)
{
}

void TopologyThreadFunction<x86_intel>::operator() ()
{
    CONTEXT("When executing TopologyThreadFunction<x86_intel> on LPU with sched / mask id " + stringify(lpu->sched_id));

    cpu_set_t _affinity_mask;
    CPU_ZERO(&_affinity_mask);
    CPU_SET(lpu->sched_id, &_affinity_mask);

    if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask) != 0)
        throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

    x86_unit * unit = retrieve_topology<x86_intel>(lpu->sched_id);
    apic_id = unit->apic_id;
    lpu->smt_id = unit->topo_id[0];
    lpu->core_id = unit->topo_id[1];
    lpu->socket_id = unit->topo_id[2];

#ifdef DEBUG
    std::string msg = "SCHED \t\t APIC \t\t SMT \t\t CORE \t\t PKG \n";
    msg += stringify(lpu->sched_id) + "\t\t" + stringify(apic_id) + "\t\t" + stringify(lpu->smt_id) + "\t\t" + stringify(lpu->core_id) +  "\t\t" + stringify(lpu->socket_id) + " \n";
    LOGMESSAGE(lc_backend, msg);
#endif

    barrier->wait();
}
#endif

template class InstantiationPolicy<Topology, Singleton>;

void Topology::determine_arch()
{
#if defined(__i386__) || defined(__x86_64__)

    /// Read processor specific information using cpuid instruction
    CONTEXT("When retrieving x86-specific vendor information from a processor via cpuid:");

    int CPUInfo[4];
    cpuid(CPUInfo, 0x0);

    if (CPUInfo[0] > 0)
    {
        if (CPUInfo[2] == 1818588270)
            _arch = x86_intel;
        else if (CPUInfo[2] == 1145913699)
            _arch = x86_amd;
        else
            _arch = unknown;
    }
#else
    _arch = unknown;
#endif
}

#if defined(__i386__) || defined(__x86_64__)
void Topology::enumerate_x86_intel()
{
    _num_cores = 0;

    Barrier * barrier = new Barrier(_num_lpus + 1);

    Thread ** threads = new Thread * [_num_lpus];

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        TopologyThreadFunction<x86_intel> ttf(_lpus[i], barrier);
        Thread * thread = new Thread(ttf);

        threads[i] = thread;
    }

    barrier->wait();
    delete barrier;

    int socket_lpu_count[_num_lpus];
    // renaming possibly necessary since socket-ids may not be continuous
    int new_socket_id[_num_lpus];

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        socket_lpu_count[i] = 0;
        new_socket_id[i] = 0;
    }

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        ++socket_lpu_count[_lpus[i]->socket_id];
        delete threads[i];
    }

    delete[] threads;

    int num_sockets(0);
    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        if (socket_lpu_count[i] > 0)
        {
            new_socket_id[i] = num_sockets;
            ++num_sockets;
        }
    }

    _num_cpus = num_sockets;

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        _lpus[i]->socket_id = new_socket_id[_lpus[i]->socket_id];
    }

    _sockets = new Socket * [num_sockets];

    int socket_lpus_assigned[num_sockets];

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        if (socket_lpu_count[i] > 0)
        {
            _sockets[new_socket_id[i]] = new Socket(socket_lpu_count[i]);
            socket_lpus_assigned[new_socket_id[i]] = 0;
        }
    }

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        int socket = _lpus[i]->socket_id;
        _sockets[socket]->_lpus[socket_lpus_assigned[socket]] = _lpus[i];
        ++socket_lpus_assigned[socket];
    }

    for (int s(0) ; s < num_sockets ; ++s)
    {
        int socket_cores[_num_lpus];

        for (unsigned l(0) ; l < _num_lpus ; ++l)
        {
            socket_cores[l] = 0;
        }

        for (int l(0) ; l < _sockets[s]->_num_lpus ; ++l)
        {
            LPU * lpu = _sockets[s]->_lpus[l];
            ++socket_cores[lpu->core_id];
        }

        int num(0);

        for (unsigned l(0) ; l < _num_lpus ; ++l)
        {
            if (socket_cores[l] > 0)
                ++num;
        }

        _sockets[s]->_num_cores = num;

        _num_cores += num;
    }

    LPU * last = _sockets[0]->_lpus[0], * lpu = _sockets[0]->_lpus[0];
    for (int l(1) ; l < _sockets[0]->_num_lpus ; ++l)
    {
       lpu = _sockets[0]->_lpus[l];
       last->linear_succ = lpu;
       last = lpu;
    }

    for (int s(1) ; s < num_sockets ; ++s)
    {
        for (int l(0) ; l < _sockets[s]->_num_lpus ; ++l)
        {
            lpu = _sockets[s]->_lpus[l];
            last->linear_succ = lpu;
            last = lpu;
        }
    }

    lpu->linear_succ = _sockets[0]->_lpus[0];

    LPU * l = _sockets[0]->_lpus[0];
    LPU * m = l;
    int pos(0);

    while (true)
    {
        int cur_sock = l->socket_id;

        if (cur_sock != num_sockets - 1)
        {
            ++cur_sock;
            m = _sockets[cur_sock]->_lpus[pos];
            l->alternating_succ = m;
        }
        else
        {
            cur_sock = 0;
            ++pos;

            if (pos < _sockets[cur_sock]->_num_lpus)
            {
                m = _sockets[cur_sock]->_lpus[pos];
                l->alternating_succ = m;
            }
            else
            {
                l->alternating_succ = _sockets[0]->_lpus[0];
                break;
            }
        }

        l = m;
    }

    int CPUInfo[4];
    cpuid(CPUInfo, 0x1);

    if ((((CPUInfo[3] >> 28) & 0x1) == 1) && (((CPUInfo[1] >> 16) & 0xFF) > 1))
        _ht_support = true;
    else
        _ht_support = false;

    if (_num_cores == _num_lpus)
        _ht_factor = 1;
    else
        _ht_factor = _num_lpus / _num_cores;

#ifdef DEBUG
    std::string msg = "Found INTEL processor(s) with " + stringify(num_sockets) + " sockets, "
        + stringify(_num_cores) + " cores and HTT " + (_ht_support == 1 ? "supported" : "not supported")
        + " and " + (_ht_factor > 1 ? "enabled" : "disabled") + "\n";

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        LPU * lpu = _lpus[i];
        msg += "LPU " + stringify(lpu->sched_id) + " on socket " + stringify(lpu->socket_id) + "\n";
        msg += "Linear succ of " + stringify(lpu->sched_id) + ": " + stringify(lpu->linear_succ->sched_id) + "\n";
        msg += "Alternating succ of " + stringify(lpu->sched_id) + ": " + stringify(lpu->alternating_succ->sched_id) + "\n";
    }

    if (_ht_factor > 1 && _ht_support == false)
        msg += "Attention: Error in SMT detection!\n";

    LOGMESSAGE(lc_backend, msg);
#endif

}

void Topology::enumerate_x86_amd()
{
    /* Topology extraction is not yet very evolved in AMD
     * processors. Try to gain information via numainfo.
     * If none is available (returned number of nodes is zero),
     * then do a simple enumeration strategy. */

    int num_nodes = intern::num_nodes();

    if (num_nodes > 0)
    {
        enumerate_numainfo(num_nodes);
    }
    else
    {
        // Assume one node - maybe done better when finding
        // the right way of comparing _num_cores with _num_lpus
        // depending on the htt/cmp_legacy values!
        _num_cpus = 1;

        int CPUInfo88[4];
        cpuid(CPUInfo88, 0x80000008);
/*
        int CPUInfo1[4];
        cpuid(CPUInfo1, 0x1);

        int htt_bit = (CPUInfo1[3] >> 28) & 0x1;
        int cmp_legacy_bit = (CPUInfo1[2] >> 1) & 0x1;
*/
        _num_cores = 1 + (CPUInfo88[2] & 0xFF);

        if (_num_cores < _num_lpus && (_num_lpus % _num_cores == 0))
        {
            _ht_support = 1;
            _ht_factor = _num_lpus / _num_cores;
        }
        else
        {
            _ht_support = 0;
            _ht_factor = 1;
        }

        _sockets = new Socket*[1];

        _sockets[0] = new Socket(_num_lpus);
        _sockets[0]->_num_cores = _num_lpus;

        for (unsigned c(0) ; c < _num_lpus ; ++c)
        {
            _sockets[0]->_lpus[c] = _lpus[c];
            _lpus[c]->socket_id = 0;
        }


        Socket * so = _sockets[0];
        LPU * last = so->_lpus[0], * lpu(0);

        for (int l(1) ; l < so->_num_lpus ; ++l)
        {
            lpu = so->_lpus[l];
            last->linear_succ = lpu;
            last->alternating_succ = lpu;
            last = lpu;
        }

        lpu->linear_succ = so->_lpus[0];
        lpu->alternating_succ = so->_lpus[0];
    }

#ifdef DEBUG
        std::string msg = "Found AMD processor(s) with " + stringify(_num_cpus) + " sockets, "
            + stringify(_num_lpus) + " LPUs and HTT " + (_ht_factor > 1 ? "available" : "not available") + "\n";

        msg += "Attention: currently no working SMT detection for AMD - always assuming #cores = #LPUs\n";

        for (unsigned i(0) ; i < _num_lpus ; ++i)
        {
            LPU * lpu = _lpus[i];
            msg += "LPU " + stringify(lpu->sched_id) + " on socket " + stringify(lpu->socket_id) + "\n";
            msg += "Linear succ of " + stringify(lpu->sched_id) + ": " + stringify(lpu->linear_succ->sched_id) + "\n";
            msg += "Alternating succ of " + stringify(lpu->sched_id) + ": " + stringify(lpu->alternating_succ->sched_id) + "\n";
        }

        LOGMESSAGE(lc_backend, msg);
#endif
}
#endif


void Topology::enumerate_numainfo(int _num_nodes)
{
    if (_num_nodes == 1)
    {
        _num_cpus = 1;
        _num_cores = _num_lpus; // Maybe possible to read SMT info from numa file system

        _sockets = new Socket*[1];

        _sockets[0] = new Socket(_num_lpus);
        _sockets[0]->_num_cores = _num_lpus;

        for (unsigned c(0) ; c < _num_lpus ; ++c)
        {
            _sockets[0]->_lpus[c] = _lpus[c];
            _lpus[c]->socket_id = 0;
        }

        Socket * so = _sockets[0];
        LPU * last = so->_lpus[0], * lpu(0);

        for (int l(1) ; l < so->_num_lpus ; ++l)
        {
            lpu = so->_lpus[l];
            last->linear_succ = lpu;
            last->alternating_succ = lpu;
            last = lpu;
        }

        lpu->linear_succ = so->_lpus[0];
        lpu->alternating_succ = so->_lpus[0];
    }
    else
    {
        unsigned * lpu_to_node = intern::cpu_to_node_array(_num_nodes, _num_lpus);

        _num_cpus = _num_nodes;
        _num_cores = _num_lpus; // Maybe possible to read SMT info from numa file system

        _sockets = new Socket*[_num_cpus];

        int socket_lpu_count[_num_lpus];

        for (unsigned i(0) ; i < _num_lpus ; ++i)
            socket_lpu_count[i] = 0;

        for (unsigned i(0) ; i < _num_lpus ; ++i)
        {
            ++socket_lpu_count[lpu_to_node[i]];
        }

        int socket_lpus_assigned[_num_cpus];

        for (unsigned i(0) ; i < _num_cpus ; ++i)
        {
            _sockets[i] = new Socket(socket_lpu_count[i]);
            socket_lpus_assigned[i] = 0;
        }

        for (unsigned i(0) ; i < _num_lpus ; ++i)
        {
            int socket = lpu_to_node[i];
            _sockets[socket]->_lpus[socket_lpus_assigned[socket]] = _lpus[i];
            _lpus[i]->socket_id = socket;
            ++socket_lpus_assigned[socket];
        }

        LPU * last = _sockets[0]->_lpus[0], * lpu(0);

        for (int l(1) ; l < _sockets[0]->_num_lpus ; ++l)
        {
           lpu = _sockets[0]->_lpus[l];
           last->linear_succ = lpu;
           last = lpu;
        }

        for (unsigned s(1) ; s < _num_cpus ; ++s)
        {
            for (int l(0) ; l < _sockets[s]->_num_lpus ; ++l)
            {
                lpu = _sockets[s]->_lpus[l];
                last->linear_succ = lpu;
                last = lpu;
            }
        }

        lpu->linear_succ = _sockets[0]->_lpus[0];

        LPU * l = _sockets[0]->_lpus[0];
        LPU * m = l;
        int pos(0);

        while (true)
        {
            unsigned cur_sock = l->socket_id;

            if (cur_sock != _num_cpus - 1)
            {
                ++cur_sock;
                m = _sockets[cur_sock]->_lpus[pos];
                l->alternating_succ = m;
            }
            else
            {
                cur_sock = 0;
                ++pos;

                if (pos < _sockets[cur_sock]->_num_lpus)
                {
                    m = _sockets[cur_sock]->_lpus[pos];
                    l->alternating_succ = m;
                }
                else
                {
                    l->alternating_succ = _sockets[0]->_lpus[0];
                    break;
                }
            }

            l = m;
        }
    }

    _ht_support = 0;
    _ht_factor = 1;

}


Topology::Topology() :
    _num_lpus(1)
{
    CONTEXT("When investigating the system topology:");

#if defined linux
    _num_lpus = sysconf(_SC_NPROCESSORS_CONF);
#endif

#ifdef DEBUG
    {
        std::string msg = "Found " + stringify(_num_lpus) + " logical processing units \n";
        LOGMESSAGE(lc_backend, msg);
    }
#endif

    _lpus = new LPU * [_num_lpus];

    for (unsigned i(0) ; i < _num_lpus ; ++i)
    {
        _lpus[i] = new LPU(i);
    }

    determine_arch();

    // Different methods must set _num_cpus, _num_cores

    switch (_arch)
    {
#if defined(__i386__) || defined(__x86_64__)
        case x86_intel:
            enumerate_x86_intel();
            break;

        case x86_amd:
            enumerate_x86_amd();
            break;
#endif
        default:
        {
            enumerate_numainfo(1);

#ifdef DEBUG
            std::string msg = "Could not determine vendor and further information about available processors\n";
            LOGMESSAGE(lc_backend, msg);
#endif
        }
    }
}


Topology::~Topology()
{
    for (unsigned c(0) ; c < _num_cpus ; ++c)
        delete _sockets[c];

    delete[] _sockets;

    for (unsigned l(0) ; l < _num_lpus ; ++l)
        delete _lpus[l];

    delete[] _lpus;
}

LPU ** Topology::lpus()
{
    return _lpus;
}

LPU * Topology::lpu(int id)
{
    return _lpus[id];
}

Socket ** Topology::sockets()
{
    return _sockets;
}

unsigned Topology::num_lpus() const
{
    return _num_lpus;
}

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
