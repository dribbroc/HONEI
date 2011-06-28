/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Sven Mallach <mallach@honei.org>
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
#ifndef HONEI_GUARD_HONEI_BACKENDS_MULTICORE_CASDEQUE_HH
#define HONEI_GUARD_HONEI_BACKENDS_MULTICORE_CASDEQUE_HH 1

#include <honei/backends/multicore/thread_task.hh>
#include <honei/util/mutex.hh>

/* Attention: This class requires the template parameter
 * to be a pointer type in order to function correctly.
 * For the moment, it's purpose is mainly to serve as a
 * CASDeque<ThreadTask *>. */

namespace honei
{
    namespace mc
    {
        namespace intern
        {
            template <typename T> struct CASDequeElement
            {
                CASDequeElement<T> * _prev;
                T _data;
                CASDequeElement<T> * _next;

                CASDequeElement() :
                    _prev(NULL),
                    _data(0),
                    _next(NULL)
                {
                }

                CASDequeElement(T & d) :
                    _prev(NULL),
                    _data(d),
                    _next(NULL)
                {
                }
            };

            template <typename T> struct CASDequeEndElement
            {
                CASDequeElement<T> * volatile _prev;
                CASDequeElement<T> * volatile _elem;
                CASDequeElement<T> * volatile _next;
                volatile int _blocked;

                CASDequeEndElement() :
                    _prev(NULL),
                    _elem(NULL),
                    _next(NULL),
                    _blocked(0)
                {
                }
            };

        }

        template <typename T> class CASDeque
        {
            private:
                intern::CASDequeEndElement<T> * const _head;
                intern::CASDequeEndElement<T> * const _tail;

            public:

                CASDeque();
                ~CASDeque();

//                void print();
                void push_back(T & data);
                T pop_front();
                T pop_back();
        };
    }
}
#endif
