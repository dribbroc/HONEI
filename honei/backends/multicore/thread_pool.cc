/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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
    num_lpus(sysconf(_SC_NPROCESSORS_CONF)),
    num_threads(Configuration::instance()->get_value("mc::num_threads", sysconf(_SC_NPROCESSORS_CONF))),
    mutex(new Mutex),
    global_barrier(new ConditionVariable),
    affinity(Configuration::instance()->get_value("mc::affinity", true))
{
    for (unsigned i(0) ; i < num_threads ; ++i)
    {
        ThreadFunction * tobj = new ThreadFunction(mutex, global_barrier, &tasks);
        Thread * t = new Thread(*tobj);
        threads.push_back(std::make_pair(t, tobj));
    }

    if (affinity)
    {
#ifdef linux
        affinity_mask = new cpu_set_t[num_threads + 1];

        // set own affinity first
        CPU_ZERO(&affinity_mask[num_threads]);
        CPU_SET(num_lpus - 1, &affinity_mask[num_threads]);
        if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &affinity_mask[num_threads]) != 0)
            throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

        unsigned i(num_threads - 1);
        for (std::list<std::pair<Thread *, ThreadFunction *> >::iterator k(threads.begin()), k_end(threads.end()) ; k != k_end ; ++k, --i)
        {
            ThreadFunction * t = (*k).second;
            thread_ids.push_back(t->tid());
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(i % (num_lpus - 1), &affinity_mask[i]);
            if(sched_setaffinity(t->tid(), sizeof(cpu_set_t), &affinity_mask[i]) != 0)
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
    if (affinity)
    {
#ifdef linux
        cpu_set_t * aff_mask = new cpu_set_t[num_threads + num + 1];
        std::copy(affinity_mask, affinity_mask + num_threads + 1, aff_mask);

        delete[] affinity_mask;
        affinity_mask = aff_mask;

        for (unsigned i(num_threads + num - 1) ; i >= num_threads ; --i)
        {
            ThreadFunction * tobj = new ThreadFunction(mutex, global_barrier, &tasks);
            Thread * t = new Thread(*tobj);
            threads.push_back(std::make_pair(t, tobj));

            thread_ids.push_back(tobj->tid());
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(i % (num_lpus - 1), &affinity_mask[i]);
            if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

        }
#endif
    }
    else
    {
        for (unsigned i(0) ; i < num ; ++i)
        {
            ThreadFunction * tobj = new ThreadFunction(mutex, global_barrier, &tasks);
            Thread * t = new Thread(*tobj);
            threads.push_back(std::make_pair(t, tobj));
        }
    }

    num_threads += num;
}

void ThreadPool::delete_threads(const unsigned num)
{
    for (unsigned i(0) ; i < num ; ++i)
    {
        std::pair<Thread *, ThreadFunction *> t = threads.back();
        threads.pop_back();
        thread_ids.pop_back();
        t.second->stop();
        delete t.second;
        delete t.first;
    }

    num_threads -= num;

    cpu_set_t * aff_mask = new cpu_set_t[num_threads + 1];
    std::copy(affinity_mask, affinity_mask + num_threads + 1, aff_mask);
    delete[] affinity_mask;
    affinity_mask = aff_mask;
}

unsigned ThreadPool::get_num_threads() const
{
    return num_threads;
}

std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > & ThreadPool::enqueue(const std::tr1::function<void ()> & task, DispatchPolicy p)
{
    std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > & ticket((affinity ? p.apply(thread_ids) : DispatchPolicy::any_core().apply(thread_ids)));

    ThreadTask * t_task(new ThreadTask(task, ticket));

    Lock l(*mutex);
    tasks.push_back(t_task);

    global_barrier->broadcast();

    return ticket;
}
