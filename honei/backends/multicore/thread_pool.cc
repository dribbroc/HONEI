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
    topology(Topology::instance()),
    num_threads(Configuration::instance()->get_value("mc::num_threads", topology->num_lpus())),
    inst_ctr(0),
    mutex(new Mutex),
    global_barrier(new ConditionVariable),
    affinity(Configuration::instance()->get_value("mc::affinity", true))
{
#ifdef linux
    if (affinity)
    {
        affinity_mask = new cpu_set_t[num_threads + 1];

        // set own affinity first
        CPU_ZERO(&affinity_mask[num_threads]);
        CPU_SET(topology->num_lpus() - 1, &affinity_mask[num_threads]);
        if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &affinity_mask[num_threads]) != 0)
            throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
    }
#endif

    for (int i(num_threads - 1) ; i >= 0 ; --i)
    {
        unsigned sched_id(i % (topology->num_lpus() - 1));

        // ToDo: Remove branch!
        ThreadFunction * tobj = new ThreadFunction(mutex, global_barrier, &tasks, inst_ctr, (affinity ? sched_id : 0xFFFF));
        Thread * t = new Thread(*tobj);
        threads.push_back(std::make_pair(t, tobj));
        ++inst_ctr;

#ifdef linux
        if (affinity)
        {
            sched_ids.push_back(sched_id);
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(sched_id, &affinity_mask[i]);
            if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
        }
#endif
    }
}

ThreadPool::~ThreadPool()
{
    for(std::list<std::pair<Thread *, ThreadFunction *> >::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        (*i).second->stop();
        delete (*i).second;
        delete (*i).first;
    }

    delete[] affinity_mask;
    delete global_barrier;
    delete mutex;
}

void ThreadPool::add_threads(const unsigned num)
{
#ifdef linux
    if (affinity)
    {
        cpu_set_t * aff_mask = new cpu_set_t[num_threads + num + 1];
        std::copy(affinity_mask, affinity_mask + num_threads + 1, aff_mask);

        delete[] affinity_mask;
        affinity_mask = aff_mask;
    }
#endif

    for (unsigned i(num_threads + num - 1) ; i >= num_threads ; --i)
    {
        unsigned sched_id(i % (topology->num_lpus() - 1));

        // ToDo: Remove branch!
        ThreadFunction * tobj = new ThreadFunction(mutex, global_barrier, &tasks, inst_ctr, (affinity ? sched_id : 0xFFFF));
        Thread * t = new Thread(*tobj);
        threads.push_back(std::make_pair(t, tobj));
        ++inst_ctr;

#ifdef linux
        if (affinity)
        {
            sched_ids.push_back(sched_id);
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(sched_id, &affinity_mask[i]);
            if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
        }
#endif
    }

    num_threads += num;
}

void ThreadPool::delete_threads(const unsigned num)
{
    for (unsigned i(0) ; i < num ; ++i)
    {
        std::pair<Thread *, ThreadFunction *> t = threads.back();
        threads.pop_back();
        sched_ids.pop_back();
        t.second->stop();
        delete t.second;
        delete t.first;
        --inst_ctr;
    }

    num_threads -= num;

    cpu_set_t * aff_mask = new cpu_set_t[num_threads + 1];
    std::copy(affinity_mask, affinity_mask + num_threads + 1, aff_mask);
    delete[] affinity_mask;
    affinity_mask = aff_mask;
}

unsigned ThreadPool::num_nodes() const
{
    return topology->num_nodes();
}

unsigned ThreadPool::get_num_threads() const
{
    return num_threads;
}

Ticket<tags::CPU::MultiCore> * ThreadPool::enqueue(const std::tr1::function<void ()> & task, DispatchPolicy p)
{
    Ticket<tags::CPU::MultiCore> * ticket((affinity ? p.apply(sched_ids) : DispatchPolicy::any_core().apply(sched_ids)));

    ThreadTask * t_task(new ThreadTask(task, ticket));

    Lock l(*mutex);
    tasks.push_back(t_task);

    global_barrier->broadcast();

    return ticket;
}
