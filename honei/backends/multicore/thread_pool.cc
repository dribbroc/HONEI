/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009, 2010 Sven Mallach <mallach@honei.org>
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
#include <honei/util/stringify.hh>

#include <errno.h>
#include <sched.h>
#include <stdlib.h>
#include <sys/syscall.h>

using namespace honei;
using namespace honei::mc;

template class InstantiationPolicy<ThreadPool, Singleton>;

ThreadPool::ThreadPool() :
    _topology(Topology::instance()),
    _num_threads(Configuration::instance()->get_value("mc::num_threads", _topology->num_lpus() - 1)),
    _inst_ctr(0),
    _mutex(new Mutex),
    _global_barrier(new ConditionVariable),
    _affinity(Configuration::instance()->get_value("mc::affinity", true))
{
#ifdef linux
    if (_affinity)
    {
        _affinity_mask = new cpu_set_t[_num_threads + 1];

        // set own affinity first
        CPU_ZERO(&_affinity_mask[_num_threads]);
        CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
        if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
            throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
    }
#endif

    for (int i(_num_threads - 1) ; i >= 0 ; --i)
    {
        unsigned sched_id(i % (_topology->num_lpus() - 1));

        ThreadFunction * tobj = new ThreadFunction(_mutex, _global_barrier, &_tasks, _inst_ctr, (_affinity ? sched_id : 0xFFFF));
        Thread * t = new Thread(*tobj);
        while (tobj->tid() == 0) ; // Wait until the thread is really setup / got cpu time for the first time

        _threads.push_back(std::make_pair(t, tobj));
        ++_inst_ctr;

#ifdef linux
        if (_affinity)
        {
            _sched_ids.push_back(sched_id);
            CPU_ZERO(&_affinity_mask[i]);
            CPU_SET(sched_id, &_affinity_mask[i]);
            if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
        }
#endif
    }
}

ThreadPool::~ThreadPool()
{
    for(std::list<std::pair<Thread *, ThreadFunction *> >::iterator i(_threads.begin()), i_end(_threads.end()) ; i != i_end ; ++i)
    {
        (*i).second->stop();
        delete (*i).second;
        delete (*i).first;
    }

    delete[] _affinity_mask;
    delete _global_barrier;
    delete _mutex;
}

void ThreadPool::add_threads(const unsigned num)
{
#ifdef linux
    if (_affinity)
    {
        cpu_set_t * aff_mask = new cpu_set_t[_num_threads + num + 1];
        std::copy(_affinity_mask, _affinity_mask + _num_threads + 1, aff_mask);

        delete[] _affinity_mask;
        _affinity_mask = aff_mask;
    }
#endif

    for (unsigned i(_num_threads + num - 1) ; i >= _num_threads ; --i)
    {
        unsigned sched_id(i % (_topology->num_lpus() - 1));

        ThreadFunction * tobj = new ThreadFunction(_mutex, _global_barrier, &_tasks, _inst_ctr, (_affinity ? sched_id : 0xFFFF));
        Thread * t = new Thread(*tobj);
        unsigned tid(0);
        do
        {
            // Wait until the thread is really setup. (got cpu time for the first time)
            tid = tobj->tid();
        }
        while (tid == 0);

        Lock l(*_mutex);

        _threads.push_back(std::make_pair(t, tobj));
        ++_inst_ctr;

#ifdef linux
        if (_affinity)
        {
            _sched_ids.push_back(sched_id);
            CPU_ZERO(&_affinity_mask[i]);
            CPU_SET(sched_id, &_affinity_mask[i]);
            if(sched_setaffinity(tid, sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
        }
#endif
    }

    _num_threads += num;
}

void ThreadPool::delete_threads(const unsigned num)
{
    for (unsigned i(0) ; i < num ; ++i)
    {
        std::pair<Thread *, ThreadFunction *> t = _threads.back();
        {
            Lock l(*_mutex);
            _sched_ids.pop_back();
            _threads.pop_back();
            t.second->stop();

        }

        delete t.second;
        delete t.first;
        --_inst_ctr;
    }

    _num_threads -= num;

#ifdef linux
    if (_affinity)
    {
        cpu_set_t * aff_mask = new cpu_set_t[_num_threads + 1];
        std::copy(_affinity_mask, _affinity_mask + _num_threads + 1, aff_mask);
        delete[] _affinity_mask;
        _affinity_mask = aff_mask;
    }
#endif
}

unsigned ThreadPool::num_nodes() const
{
    return _topology->num_nodes();
}

unsigned ThreadPool::num_threads() const
{
    return _num_threads;
}

unsigned ThreadPool::main_node() const
{
    return _topology->get_node(_topology->num_lpus() - 1);
}

Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const function<void ()> & task, DispatchPolicy p)
{
    Ticket<tags::CPU::MultiCore> * ticket((_affinity ? p.apply(_sched_ids) : DispatchPolicy::any_core().apply(_sched_ids)));

    ThreadTask * t_task(new ThreadTask(task, ticket));

    Lock l(*_mutex);
    _tasks.push_back(t_task);

    _global_barrier->broadcast();

    return ticket;
}
