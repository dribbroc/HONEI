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
#include <honei/util/configuration.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/stringify.hh>

//todo wenn das alles fertig ist, scaledsum auf multigpu

#include <errno.h>
#include <sched.h>
#include <stdlib.h>
#include <sys/syscall.h>

using namespace honei;
using namespace honei::cuda;

template class InstantiationPolicy<GPUPool, Singleton>;

GPUPool::GPUPool() :
    num_gpus(2),
    mutex(new Mutex),
    global_barrier(new ConditionVariable)
{
    //for (unsigned i(0) ; i < num_gpus ; ++i)
    //{
    // todo auf task vector tasks[device] umstellen
        GPUFunction * tobj = new GPUFunction(0, mutex, global_barrier, &tasks_gpu0);
        Thread * t = new Thread(*tobj);
        threads.push_back(std::make_pair(t, tobj));
        GPUFunction * tobj2 = new GPUFunction(1, mutex, global_barrier, &tasks_gpu1);
        Thread * t2 = new Thread(*tobj2);
        threads.push_back(std::make_pair(t2, tobj2));
    //}
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
    // todo auf task vector tasks[device] umstellen und exception wenn falsches device angegeben
    switch (device)
    {
        case 0:
            tasks_gpu0.push(t_task);
            break;
        case 1:
            tasks_gpu1.push(t_task);
            break;
    }


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
    // todo auch checken ob die queues leer sind
    return true;
}

void GPUPool::flush()
{
    SynchTask t;
    TicketVector tickets;
    tickets.push_back(enqueue(t,0));
    tickets.push_back(enqueue(t,1));
    tickets.wait();
}
