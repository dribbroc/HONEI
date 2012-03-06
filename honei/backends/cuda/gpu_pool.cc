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
    num_gpus(std::min(Configuration::instance()->get_value("cuda::num_gpu", 2), cuda_device_count()))
{
    if (num_gpus == 0)
        throw InternalError("No GPU Found!");
    for (unsigned i(0) ; i < num_gpus ; ++i)
    {
        tasks.push_back(new std::queue<GPUTask *>);
        barriers.push_back(new ConditionVariable);
        mutexe.push_back(new Mutex);
        GPUFunction * tobj(0);
        tobj = new GPUFunction(i, mutexe.at(i), barriers.at(i), (tasks.at(i)));
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

    for (unsigned i(0) ; i < mutexe.size() ; ++i)
    {
        delete barriers.at(i);
        delete mutexe.at(i);
        delete tasks.at(i);
    }
}

unsigned GPUPool::get_num_gpus() const
{
    return num_gpus;
}

Ticket<tags::GPU::MultiCore> GPUPool::enqueue(const function<void ()> & task, int device)
{
    Ticket<tags::GPU::MultiCore> ticket;

    GPUTask * t_task(new GPUTask(task, ticket));

    Lock l(*mutexe.at(device%num_gpus));
    tasks.at(device%num_gpus)->push(t_task);

    barriers.at(device%num_gpus)->broadcast();

    return ticket;
}

bool GPUPool::idle()
{
    for(std::vector<std::pair<Thread *, GPUFunction *> >::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        if (!(*i).second->idle())
            return false;
    }
    for (unsigned i(0) ; i < num_gpus ; ++i)
    {
        Lock l(*mutexe.at(i));
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

void GPUPool::single_start(int device)
{
    if (!idle())
        throw InternalError("Cannot restart in single mode if GPUPool is not idle!");

    if (device >= cuda_device_count())
        throw InternalError("Invalid device!");

    for(std::vector<std::pair<Thread *, GPUFunction *> >::iterator i(threads.begin()), i_end(threads.end()) ; i != i_end ; ++i)
    {
        (*i).second->stop();
        delete (*i).second;
        delete (*i).first;
    }

    threads.clear();
    num_gpus = 1;

    GPUFunction * tobj = new GPUFunction(device, mutexe.at(0), barriers.at(0), (tasks.at(0)));
    Thread * t = new Thread(*tobj);
    threads.push_back(std::make_pair(t, tobj));

    cuda_set_device(device);
}
