/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/stringify.hh>

#include <errno.h>
#include <sched.h>
#include <stdlib.h>

using namespace honei;
using namespace honei::mc;

ThreadPool::ThreadPool(const unsigned n_threads, const bool aff) :
    num_lpus(sysconf(_SC_NPROCESSORS_CONF)),
    num_threads(n_threads),
    mutex(new Mutex),
    global_barrier(new ConditionVariable),
    thread_ids(new unsigned[n_threads]),
    affinity(aff)
{
    for (unsigned i(0) ; i < num_threads ; ++i)
        threads.push_back(new Thread(mutex, &tasks, global_barrier));

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
        for (std::list<Thread *>::iterator k(threads.begin()), k_end(threads.end()) ; k != k_end ; ++k, --i)
        {
            Thread * t = *k;
            thread_ids[i] = t->tid();
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(i % (num_lpus - 1), &affinity_mask[i]);
            if(sched_setaffinity(thread_ids[i], sizeof(cpu_set_t), &affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

        }
#endif
    }
}

ThreadPool * ThreadPool::instance(unsigned threads, bool affinity)
{
    if (ThreadPool::instance_ptr == 0)
    {
        void * address(0);
        posix_memalign(&address, 128, sizeof(ThreadPool));
        ThreadPool::instance_ptr = new(address) ThreadPool(threads, affinity);
    }

    return ThreadPool::instance_ptr;
}

void ThreadPool::destroy()
{
    ThreadPool::instance_ptr->~ThreadPool();
    free(ThreadPool::instance_ptr);
    ThreadPool::instance_ptr = 0;
}

ThreadPool::~ThreadPool()
{
    for(std::list<Thread *>::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        delete *i;
    }

    delete[] thread_ids;
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

        unsigned * ids = new unsigned[num_threads + num];
        std::copy(thread_ids, thread_ids + num_threads, ids);
        delete[] thread_ids;
        thread_ids = ids;

        for (unsigned i(num_threads + num - 1) ; i >= num_threads ; --i)
        {
            Thread * t = new Thread(mutex, &tasks, global_barrier);
            threads.push_back(t);
            thread_ids[i] = t->tid();
            CPU_ZERO(&affinity_mask[i]);
            CPU_SET(i % (num_lpus - 1), &affinity_mask[i]);
            if(sched_setaffinity(thread_ids[i], sizeof(cpu_set_t), &affinity_mask[i]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

        }
#endif
    }
    else
    {
        for (unsigned i(0) ; i < num ; ++i)
            threads.push_back(new Thread(mutex, &tasks, global_barrier));
    }

    num_threads += num;
}


// ToDo: We have a problem when s.o. deletes threads here that he later wants to reuse per same_core_as()
void ThreadPool::delete_threads(const unsigned num)
{
    for (unsigned i(0) ; i < num ; ++i)
    {
        Thread * t = threads.back();
        threads.pop_back();
        delete t;
    }

    num_threads -= num;

    cpu_set_t * aff_mask = new cpu_set_t[num_threads + 1];
    std::copy(affinity_mask, affinity_mask + num_threads + 1, aff_mask);
    delete[] affinity_mask;
    affinity_mask = aff_mask;

    unsigned * ids = new unsigned[num_threads];
    std::copy(thread_ids, thread_ids + num_threads, ids);
    delete[] thread_ids;
    thread_ids = ids;
}

unsigned ThreadPool::get_num_threads() const
{
    return num_threads;
}

Ticket & ThreadPool::enqueue(const std::tr1::function<void ()> & task, DispatchPolicy p)
{
    Ticket * ticket = new Ticket;

    unsigned & thread_id = p.apply(ticket, thread_ids);

    ThreadTask * t_task = new ThreadTask(task, ticket, &thread_id);

    Lock l(*mutex);
    tasks.push_back(t_task);

    global_barrier->broadcast();

    return *ticket;
}

ThreadPool * ThreadPool::instance_ptr = 0;
