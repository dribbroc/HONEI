/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012 Sven Mallach <mallach@honei.org>
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

#pragma once
#ifndef HONEI_GUARD_HONEI_BACKENDS_MULTICORE_LPU_HH
#define HONEI_GUARD_HONEI_BACKENDS_MULTICORE_LPU_HH 1

struct LPU
{
    // Will always be set
    const int sched_id; // scheduler / cpuinfo / mask id
    int socket_id; // socket / cpu / node id

    // Will not always be set
    int smt_id; // hardware thread id
    int core_id; // core id

    LPU * linear_succ;
    LPU * alternating_succ;

    // Information that will be set by the thread pool
    bool has_thread;
    LPU * linear_enqueue_succ;
    LPU * alternating_enqueue_succ;

    LPU(const int sid) :
        sched_id(sid),
        socket_id(-1),
        smt_id(-1),
        core_id(-1),
        has_thread(false)
    {
    }
};

struct Socket
{
    int _num_lpus;
    int _num_cores;
    LPU ** _lpus;
    int _num_threads;
    LPU ** _threaded_lpus;

    Socket(int num_lpus) :
        _num_lpus(num_lpus),
        _num_cores(0),
        _lpus(new LPU*[num_lpus]),
        _num_threads(0),
        _threaded_lpus(new LPU*[num_lpus])
    {
    }

    ~Socket()
    {
        delete[] _lpus;
        delete[] _threaded_lpus;
    }
};

#endif
