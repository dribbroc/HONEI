/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009, 2010, 2011 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/thread_pool.hh>
#include <honei/util/configuration.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>

#include <errno.h>
#include <sched.h>
#include <stdlib.h>
#include <sys/syscall.h>

namespace honei
{
    template <> struct Implementation<mc::ThreadPool>
    {
        /// Information about processor topology (such as number of processing units)
        mc::Topology * const _topology;

        /// Number of demanded threads (via honeirc)
        const unsigned _demanded_threads;

        /// Number of pooled threads
        const unsigned _num_threads;

        /// A thread instantiation counter
        unsigned _inst_ctr;

        /// List of user POSIX threads
        std::list<std::pair<Thread *, mc::ThreadFunctionBase *> > _threads;

        /// Waiting list of worker tasks to be executed (if affinity is enabled)
        std::list<mc::ThreadTask *> _tasks;

        /// Waiting list of worker tasks to be executed (otherwise)
        mc::AtomicSList<mc::ThreadTask *> _ttasks;

        /// Exchange data between the pool and its threads
        mc::PoolSyncData * const _pool_sync;

        /// Flag whether to use thread affinity
        const bool _affinity;

        /// Function pointer to any of the default dispatch strategies
        mc::DispatchPolicy (* policy) ();

        /// Mapping of threads to the scheduler ids of the cores they run on
        std::vector<unsigned> _sched_ids;

#ifdef linux
        /// Array of affinity masks for main process and all controlled threads
        cpu_set_t * _affinity_mask;
#endif

        Implementation() :
            _topology(mc::Topology::instance()),
            _demanded_threads(Configuration::instance()->get_value("mc::num_threads", _topology->num_lpus())),
            _num_threads(_demanded_threads > _topology->num_lpus() ? _demanded_threads : _topology->num_lpus()),
            _inst_ctr(0),
            _pool_sync(new mc::PoolSyncData),
            _affinity(Configuration::instance()->get_value("mc::affinity", true))
        {
            CONTEXT("When initializing the thread pool:\n");

#ifndef linux
            _affinity = false;
#endif

#ifdef DEBUG
            std::string msg = "Will create " + stringify(_num_threads) + " POSIX worker threads\n";
            if (_affinity)
                msg += "Found Thread Affinity to be enabled\n";
            else
                msg += "Found Thread Affinity to be disabled\n";

            msg += "Default dispatch policy is configured as ";
#endif

            std::string dis = Configuration::instance()->get_value("mc::dispatch", "anycore");

            if (! _affinity || dis == "anycore")
            {
                policy = &mc::DispatchPolicy::any_core;
#ifdef DEBUG
                msg += "arbitrary - the next available thread will execute a task\n";
#endif
            }
            else if (dis == "alternating")
            {
                policy = &mc::DispatchPolicy::alternating_node;
#ifdef DEBUG
                msg += "alternating - tasks shall be assigned alternatingly on available NUMA nodes\n";
#endif
            }
            else if (dis == "linear")
            {
                policy = &mc::DispatchPolicy::linear_node;
#ifdef DEBUG
                msg += "linear - tasks shall be assigned to fill up available NUMA nodes one-by-one\n";
#endif
            }

            if (_affinity)
            {
                _affinity_mask = new cpu_set_t[_num_threads + 1];

                // set own affinity first
                CPU_ZERO(&_affinity_mask[_num_threads]);
                CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
                if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
                    throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
                msg += "THREAD \t\t POOL_ID \t LPU \t NODE \n";
                msg += "MAIN \t\t - \t\t" + stringify(_topology->num_lpus() - 1) + "\t\t" + stringify(_topology->get_node(_topology->num_lpus() - 1)) + " \n";
#endif
            }

            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                if (_affinity)
                {
                    unsigned sched_id(i % (_topology->num_lpus()));
                    mc::AffinityThreadFunction * tobj = new mc::AffinityThreadFunction(_pool_sync, &_tasks, _inst_ctr, sched_id);
                    Thread * t = new Thread(*tobj);
                    while (tobj->tid() == 0) ; // Wait until the thread is really setup / got cpu time for the first time
                    _threads.push_back(std::make_pair(t, tobj));
                    _sched_ids.push_back(sched_id);
                    CPU_ZERO(&_affinity_mask[i]);
                    CPU_SET(sched_id, &_affinity_mask[i]);
                    if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                        throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
#ifdef DEBUG
                    msg += stringify(tobj->tid()) + "\t\t" + stringify(_inst_ctr) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif
                }
                else
                {
                    mc::SimpleThreadFunction * tobj = new mc::SimpleThreadFunction(_pool_sync, &_ttasks, _inst_ctr);
                    Thread * t = new Thread(*tobj);
                    while (tobj->tid() == 0) ; // Wait until the thread is really setup / got cpu time for the first time
                    _threads.push_back(std::make_pair(t, tobj));
                }

                ++_inst_ctr;
            }

#ifdef DEBUG
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        ~Implementation()
        {
            for(std::list<std::pair<Thread *, mc::ThreadFunctionBase *> >::iterator i(_threads.begin()),
                i_end(_threads.end()) ; i != i_end ; ++i)
            {
                (*i).second->stop();
                delete (*i).second;
                delete (*i).first;
            }

            if (_affinity)
                delete[] _affinity_mask;

            delete _pool_sync;
        }
    };
}

using namespace honei;
using namespace honei::mc;

template class InstantiationPolicy<ThreadPool, Singleton>;

Ticket<tags::CPU::MultiCore> * DispatchPolicy::last(NULL);


ThreadPool::ThreadPool() :
    PrivateImplementationPattern<ThreadPool, Single>(new Implementation<ThreadPool>)
{
}

ThreadPool::~ThreadPool()
{
}

unsigned ThreadPool::num_threads() const
{
    return _imp->_num_threads;
}

Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const function<void ()> & task, DispatchPolicy p)
{
    CONTEXT("When creating a ThreadTask:\n");

    Ticket<tags::CPU::MultiCore> * ticket(NULL);

    if (_imp->_affinity)
    {
        ticket = p.apply();
        ThreadTask * t_task(new ThreadTask(task, ticket));
        Lock l(*_imp->_pool_sync->mutex);
        _imp->_tasks.push_back(t_task);
        _imp->_pool_sync->barrier->broadcast();
    }
    else
    {
        ticket = DispatchPolicy::any_core().apply();
        ThreadTask * t_task(new ThreadTask(task, ticket));
        _imp->_ttasks.push_back(t_task);
        Lock l(*_imp->_pool_sync->mutex);
        _imp->_pool_sync->barrier->broadcast();
    }

    return ticket;
}

/// Use default policy
Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const function<void ()> & task)
{
    CONTEXT("When creating a ThreadTask:\n");

    Ticket<tags::CPU::MultiCore> * ticket(_imp->policy().apply());

    ThreadTask * t_task(new ThreadTask(task, ticket));

    if (_imp->_affinity)
    {
        Lock l(*_imp->_pool_sync->mutex);
        _imp->_tasks.push_back(t_task);
        _imp->_pool_sync->barrier->broadcast();
    }
    else
    {
        _imp->_ttasks.push_back(t_task);
        Lock l(*_imp->_pool_sync->mutex);
        _imp->_pool_sync->barrier->broadcast();
    }

    return ticket;
}
