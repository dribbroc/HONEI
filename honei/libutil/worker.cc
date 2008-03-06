/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/libutil/exception.hh>
#include <honei/libutil/instantiation_policy-impl.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/mutex.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/worker.hh>

#include <list>

#include <cstring>
#include <cerrno>

using namespace honei;

struct WorkerThread::Implementation
{
    typedef std::list<WorkerTask> TaskList;

    /// Our thread.
    pthread_t * const thread;

    /// Our mutex for read/write access to a WorkerThread instance.
    Mutex * const mutex;

    /// Our exit condition.
    bool exit;

    /// Our list of queued tasks.
    TaskList task_list;

    /// Get the next task from the queue.
    inline WorkerTask * get_next()
    {
        Lock l(*mutex);
        WorkerTask * result(0);

        if (! task_list.empty())
            result = new WorkerTask(task_list.front());

        return result;
    }

    /// Dequeue the next task.
    inline void dequeue()
    {
        Lock l(*mutex);

        if (! task_list.empty())
            task_list.pop_front();
    }

    /// Enqueue a task with us.
    inline void enqueue(WorkerTask & task)
    {
        Lock l(*mutex);

        task_list.push_back(task);
    }

    /// Return true if we are idling.
    inline bool idle() const
    {
        Lock l(*mutex);

        return task_list.empty();
    }

    /// Our thread's main function.
    static void * thread_function(void * argument)
    {
        Implementation * imp(static_cast<Implementation *>(argument));
        WorkerTask * task(0);
        bool exit(false);

        do
        {
            // Get a task
            task = imp->get_next();

            if (! task)
            {
                sleep(1);
            }

            // Update our exit condition
            {
                Lock l(*imp->mutex);
                exit = imp->exit;
            }

            // Run our task
            if (task)
            {
                (*task)();
                imp->dequeue();
            }
        }
        while (! exit);

        pthread_exit(0);
    }

    Implementation() :
        thread(new pthread_t),
        mutex(new Mutex),
        exit(false)
    {
        int retval;

        if (0 != (retval = pthread_create(thread, 0, &thread_function, this)))
            throw ExternalError("libpthread", "pthread_create failed, " + stringify(strerror(retval)));
    }

    ~Implementation()
    {
        // Flag for exit.
        {
            Lock l(*mutex);
            exit = true;
        }

        pthread_join(*thread, 0);

        delete mutex;
        delete thread;
    }
};

WorkerThread::WorkerThread() :
    _imp(new Implementation)
{
}

WorkerThread::~WorkerThread()
{
    delete _imp;
}

void
WorkerThread::enqueue(WorkerTask & task)
{
    _imp->enqueue(task);
}

bool
WorkerThread::idle() const
{
    return _imp->idle();
}
