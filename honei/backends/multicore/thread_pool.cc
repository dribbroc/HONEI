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
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>

#include <errno.h>
#include <sched.h>
#include <stdlib.h>
#include <sys/syscall.h>

using namespace honei;
using namespace honei::mc;

template class InstantiationPolicy<ThreadPool, Singleton>;

Ticket<tags::CPU::MultiCore> * DispatchPolicy::last(NULL);

ThreadPool::ThreadPool() :
    _topology(Topology::instance()),
    _demanded_threads(Configuration::instance()->get_value("mc::num_threads", _topology->num_lpus())),
    _num_threads(_demanded_threads > _topology->num_lpus() ? _demanded_threads : _topology->num_lpus()),
    _inst_ctr(0),
    _pool_sync(new PoolSyncData),
    _affinity(Configuration::instance()->get_value("mc::affinity", true))
{
    CONTEXT("When initializing the thread pool:\n");

#ifdef DEBUG
    std::string msg = "Will create " + stringify(_num_threads) + " POSIX worker threads \n";
#endif

#ifndef linux
    _affinity = false;
#endif

    std::string dis = Configuration::instance()->get_value("mc::dispatch", "anycore");

    if (! _affinity || dis == "anycore")
        policy = &DispatchPolicy::any_core;
    else if (dis == "alternating")
        policy = &DispatchPolicy::alternating_node;
    else if (dis == "linear")
        policy = &DispatchPolicy::linear_node;

    if (_affinity)
    {
#ifdef DEBUG
        msg += "Thread Affinity is enabled - assigninh threads to definite logical processing units. \n";
#endif
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
            AffinityThreadFunction * tobj = new AffinityThreadFunction(_pool_sync, &_tasks, _inst_ctr, sched_id);
            Thread * t = new Thread(*tobj);
            while (tobj->tid() == 0) ; // Wait until the thread is really setup / got cpu time for the first time
            _threads.push_back(std::make_pair(t, tobj));
            _sched_ids.push_back(sched_id);
            CPU_ZERO(&_affinity_mask[i]);
            CPU_SET(sched_id, &_affinity_mask[i]);
            if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
            msg += stringify(tobj->tid()) + "\t\t" + stringify(_inst_ctr - 1) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif
        }
        else
        {
            SimpleThreadFunction * tobj = new SimpleThreadFunction(_pool_sync, &_ttasks, _inst_ctr);
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

ThreadPool::~ThreadPool()
{
    for(std::list<std::pair<Thread *, ThreadFunctionBase *> >::iterator i(_threads.begin()),
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

unsigned ThreadPool::num_threads() const
{
    return _num_threads;
}

Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const function<void ()> & task, DispatchPolicy p)
{
    CONTEXT("When creating a ThreadTask:\n");

    Ticket<tags::CPU::MultiCore> * ticket((_affinity ? p.apply(_sched_ids) : DispatchPolicy::any_core().apply(_sched_ids)));

    ThreadTask * t_task(new ThreadTask(task, ticket));

    if (_affinity)
    {
        Lock l(*_pool_sync->mutex);
        _tasks.push_back(t_task);
    }
    else
        _ttasks.push_back(t_task);

    _pool_sync->barrier->broadcast();

    return ticket;
}

/// Use default policy
Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const function<void ()> & task)
{
    CONTEXT("When creating a ThreadTask:\n");

    Ticket<tags::CPU::MultiCore> * ticket(policy().apply(_sched_ids));

    ThreadTask * t_task(new ThreadTask(task, ticket));

    if (_affinity)
    {
        Lock l(*_pool_sync->mutex);
        _tasks.push_back(t_task);
        _pool_sync->barrier->broadcast();
    }
    else
    {
        _ttasks.push_back(t_task);
        Lock l(*_pool_sync->mutex);
        _pool_sync->barrier->broadcast();
    }

    return ticket;
}
