/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_POOL_THREAD_HH
#define LIBUTIL_GUARD_POOL_THREAD_HH 1

#include <honei/libutil/pool_task.hh>

#include <tr1/functional>

namespace honei
{
    typedef std::tr1::function<void () throw ()> WorkerTask;

    class ThreadPool;

    class PoolThread
    {
        private:
            struct Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// \name Unwanted operations.
            /// \{

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            PoolThread(const PoolThread &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            PoolThread & operator= (const PoolThread &);

            /// \}

        public:
            /// \name Constructor and destructor
            /// \{

            /// Constructor.
            PoolThread(ThreadPool * pool);

            /// Destructor.
            ~PoolThread();

            /// \}

            /// Enqueue a WorkerTask with us.
            void run(PoolTask * task);
   };
}

#endif
