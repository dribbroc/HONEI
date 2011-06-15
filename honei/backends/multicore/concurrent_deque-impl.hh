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

#include <honei/backends/multicore/concurrent_deque.hh>

using namespace honei;
using namespace honei::mc;
using namespace honei::mc::intern;

template <typename T>
ConcurrentDeque<T>::ConcurrentDeque() :
    _head(NULL),
    _tail(NULL),
    _front(new Mutex),
    _back(new Mutex)
{
}

template <typename T>
ConcurrentDeque<T>::~ConcurrentDeque()
{
    delete _back;
    delete _front;
}


template <typename T>
bool ConcurrentDeque<T>::empty()
{
    pthread_mutex_lock(_front->mutex());
    return _head == NULL;
    pthread_mutex_unlock(_front->mutex());
}

template <typename T>
int ConcurrentDeque<T>::size()
{
    int s(0);
    DequeElement<T> * c(NULL);

    pthread_mutex_lock(_front->mutex());
    pthread_mutex_lock(_back->mutex());
    c = _head;
    while (c != NULL)
    {
        c = c->_next;
        ++s;
    }
    pthread_mutex_unlock(_back->mutex());
    pthread_mutex_unlock(_front->mutex());

    return s;
}
template <typename T>
void ConcurrentDeque<T>::push_back(T & data)
{
    DequeElement<T> * e = new DequeElement<T>(data);

    pthread_mutex_lock(_front->mutex());
    pthread_mutex_lock(_back->mutex());

    if (_head == NULL) // zero
    {
        _head = e;
        _tail = e;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        return;
    }

    if (_head->_next == NULL) // one
    {
        e->_prev = _head;
        _head->_next = e;
        _tail = e;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        return;
    }

    pthread_mutex_unlock(_front->mutex());

    // no difference between two and more

    e->_prev = _tail;
    _tail->_next = e;
    _tail = e;

    pthread_mutex_unlock(_back->mutex());
}

template <typename T>
T ConcurrentDeque<T>::pop_front()
{
    T retval(0);
    DequeElement<T> * fr(NULL);

    pthread_mutex_lock(_front->mutex());
    pthread_mutex_lock(_back->mutex());

    if (_head == NULL) // zero
    {
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        return retval;
    }

    if (_head->_next == NULL) // one
    {
        fr = _head;
        _head = NULL;
        _tail = NULL;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        retval = fr->_data;
        delete fr;
        return retval;
    }

    /* if tail = 2 becomes head then _tail->prev
     * will never be needed - so we can make no
     * difference between 2 and more elements.

    if (_head->_next == _tail) // two
    {
        fr = _head;
        _head = _head->_next;
        _tail->_prev = NULL;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        retval = fr->_data;
        delete fr;
        return retval;
    }
    */

    pthread_mutex_unlock(_back->mutex());

    fr = _head;
    _head = _head->_next;
    // no need to modify a prev pointer
    pthread_mutex_unlock(_front->mutex());
    retval = fr->_data;
    delete fr;
    return retval;
}

template <typename T>
T ConcurrentDeque<T>::pop_back()
{
    T retval(0);
    DequeElement<T> * fr(NULL);

    pthread_mutex_lock(_front->mutex());
    pthread_mutex_lock(_back->mutex());

    if (_head == NULL) // zero
    {
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        return retval;
    }

    if (_head->_next == NULL) // one
    {
        fr = _head;
        _head = NULL;
        _tail = NULL;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        retval = fr->_data;
        delete fr;
        return retval;
    }

    if (_head->_next == _tail) // two
    {
        fr = _tail;
        _head->_next = NULL;
        _tail = _head;
        pthread_mutex_unlock(_back->mutex());
        pthread_mutex_unlock(_front->mutex());
        retval = fr->_data;
        delete fr;
        return retval;
    }

    pthread_mutex_unlock(_front->mutex());

    fr = _tail;
    _tail = _tail->_prev;
    _tail->_next = NULL;
    // no need to modify a prev pointer
    pthread_mutex_unlock(_back->mutex());
    retval = fr->_data;
    delete fr;
    return retval;
}
