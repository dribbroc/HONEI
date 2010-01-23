/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@tu-dortmund.de>
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

#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/backends/cuda/multi_gpu.hh>
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
using namespace honei::cuda;

template class InstantiationPolicy<GPUPool, Singleton>;

GPUPool::GPUPool() :
    num_gpus(std::min(Configuration::instance()->get_value("cuda::num_gpu", 2), cuda_device_count())),
    mutex(new Mutex),
    global_barrier(new ConditionVariable)
{
    for (unsigned i(0) ; i < num_gpus ; ++i)
    {
        std::queue<GPUTask *> * q = new std::queue<GPUTask *>;
        tasks.push_back(q);
        GPUFunction * tobj = new GPUFunction(i, mutex, global_barrier, (tasks.at(i)));
        Thread * t = new Thread(*tobj);
        threads.push_back(std::make_pair(t, tobj));
    }
}

GPUPool::~GPUPool()
{
    for(std::vector<std::pair<Thread *, GPUFunction *> >::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        (*i).second->stop();
        delete (*i).second;
        delete (*i).first;
    }

    delete global_barrier;
    delete mutex;
}

unsigned GPUPool::get_num_gpus() const
{
    return num_gpus;
}

Ticket<tags::GPU::MultiCore> * GPUPool::enqueue(const std::tr1::function<void ()> & task, int device)
{
    Ticket<tags::GPU::MultiCore> * ticket = new Ticket<tags::GPU::MultiCore>();

    GPUTask * t_task(new GPUTask(task, ticket));

    Lock l(*mutex);
    tasks.at(device%num_gpus)->push(t_task);

    global_barrier->broadcast();

    return ticket;
}

bool GPUPool::idle()
{
    Lock l(*mutex);
    for(std::vector<std::pair<Thread *, GPUFunction *> >::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        if (!(*i).second->idle())
            return false;
    }
    for (unsigned i(0) ; i < num_gpus ; ++i)
    {
        if (tasks.at(i)->size() != 0)
            return false;
    }
    return true;
}

void GPUPool::flush()
{
    SynchTask t;
    TicketVector tickets;
    for (unsigned i(0) ; i < num_gpus ; ++i)
    {
        tickets.push_back(enqueue(t,i));
    }
    tickets.wait();
    while (!idle())
    {
    }
}
