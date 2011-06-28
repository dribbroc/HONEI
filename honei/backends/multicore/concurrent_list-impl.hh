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
#ifndef HONEI_GUARD_HONEI_BACKENDS_MULTICORE_CONCURRENT_LIST_IMPL_HH
#define HONEI_GUARD_HONEI_BACKENDS_MULTICORE_CONCURRENT_LIST_IMPL_HH 1

#include <honei/backends/multicore/concurrent_list.hh>

using namespace honei;
using namespace honei::mc;
using namespace honei::mc::intern;

template <typename T>
ConcurrentList<T>::ConcurrentList() :
    _head(NULL),
    _tail(NULL),
    _global(new Mutex),
    _front(new Mutex),
    _back(new Mutex),
    _size(0)
{
}

template <typename T>
ConcurrentList<T>::~ConcurrentList()
{
    delete _back;
    delete _front;
    delete _global;
}

template <typename T>
bool ConcurrentList<T>::acquire_front()
{
    bool ret(false);
    pthread_mutex_lock(_global->mutex());
    pthread_mutex_lock(_front->mutex());

    if (_size == 0 || _size-- < 3)
    {
        pthread_mutex_lock(_back->mutex());
        ret = true;
    }

    pthread_mutex_unlock(_global->mutex());

    return ret;
}

template <typename T>
void ConcurrentList<T>::release_front(bool back)
{
    if (back)
        pthread_mutex_unlock(_back->mutex());

    pthread_mutex_unlock(_front->mutex());
}

template <typename T>
bool ConcurrentList<T>::acquire_back()
{
    bool ret(false);
    pthread_mutex_lock(_global->mutex());
    pthread_mutex_lock(_back->mutex());

    if (_size++ < 3)
    {
        pthread_mutex_lock(_front->mutex());
        ret = true;
    }

    pthread_mutex_unlock(_global->mutex());

    return ret;
}

template <typename T>
void ConcurrentList<T>::release_back(bool front)
{
    if (front)
    {
        pthread_mutex_unlock(_front->mutex());
    }

    pthread_mutex_unlock(_back->mutex());
}

template <typename T>
bool ConcurrentList<T>::empty()
{
    pthread_mutex_lock(_global->mutex());
    bool ret = (_size == 0);
    pthread_mutex_unlock(_global->mutex());

    return ret;
}

template <typename T>
void ConcurrentList<T>::push_back(T & data)
{
    CListElement<T> * e = new CListElement<T>(data);

    bool fr = acquire_back();

    if (_tail == NULL) // Empty list;
    {
        _head = e;
        _tail = e;
    }
    else // tail and head cannot be NULL in this clause
    {
        if (fr && _tail == _head) // Only one element.
            _head->_next = e;
        else
            _tail->_next = e;

        _tail = e;
    }

    release_back(fr);
}

template <typename T>
T ConcurrentList<T>::pop_front()
{
    T retval(0);
    CListElement<T> * fr(NULL);

    bool ba = acquire_front();

    if (_head != NULL)
    {
        fr = _head;
        retval = (_head->_data);

        if (ba && _tail == _head) // Only one element
        {
            _head = NULL;
            _tail = NULL;
        }
        else
            _head = _head->_next;
    }

    release_front(ba);

    delete fr;

    return retval;
}
#endif
