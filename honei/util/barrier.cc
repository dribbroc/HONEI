/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
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

#include <honei/util/barrier.hh>

using namespace honei;

Barrier::Barrier(unsigned thread_count) :
    _barrier(new pthread_barrier_t),
    _restriction(new pthread_barrierattr_t)
{
    pthread_barrierattr_init(_restriction);
    pthread_barrier_init(_barrier, _restriction, thread_count);
}

Barrier::~Barrier()
{
    pthread_barrier_destroy(_barrier);
    pthread_barrierattr_destroy(_restriction);
    delete _barrier;
    delete _restriction;
}

void Barrier::wait()
{
    pthread_barrier_wait(_barrier);
}
