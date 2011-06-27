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

#include <honei/backends/multicore/cas_deque.hh>
#include <iostream>

using namespace honei;
using namespace honei::mc;
using namespace honei::mc::intern;

inline bool CAS(volatile int * ptr, int cmp, int val) __attribute__((always_inline));

bool CAS(volatile int * ptr, int cmp, int val)
{
    return __sync_bool_compare_and_swap(ptr, cmp, val);
}

template <typename T>
CASDeque<T>::CASDeque() :
    _head(new CASDequeEndElement<T>),
    _tail(new CASDequeEndElement<T>)
{
    std::cout << "Created list" << this << " with _head = " << _head << " and _tail = " << _tail << std::endl;
}

template <typename T>
CASDeque<T>::~CASDeque()
{
    delete _tail;
    delete _head;
}

template <typename T>
void CASDeque<T>::push_back(T & data)
{
    bool ok(false);
    CASDequeElement<T> * e = new CASDequeElement<T>(data);

    do
    {
        if (_head->_elem == NULL && _tail->_elem == NULL) // assume zero
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem == NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == NULL)
                {
                    /*
                    if (_head->_next != NULL || _tail->_prev != NULL)
                    {
                        std::cout << this << " PB: CORRUPTION ON EMPTY." << std::endl;
                    }
                    */

                    _head->_elem = e;
                    _head->_next = NULL;
                    _tail->_prev = NULL;
                    _tail->_elem = e;
//                    std::cout << this << " push-back on empty" << std::endl;
//                    print();
                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
        else if (_head->_next == NULL && _tail->_prev == NULL) // assume one
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev == NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == _tail->_elem && _head->_next == NULL)
                {
                    /*
                    if (_head->_elem == NULL || _tail->_elem == NULL)
                    {
                        std::cout << this << " PB: CORRUPTION ON ONE." << std::endl;
                    }
                    */
                    e->_prev = _tail->_elem;

                    _head->_elem->_next = e;
                    _head->_next = e;

                    _tail->_prev = _head->_elem;
                    _tail->_elem = e;
//                    std::cout << this << " push-back on one" << std::endl;
//                    print();
                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
        else if (_tail->_prev != NULL && _head->_next != NULL)
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem != NULL && _head->_elem != _tail->_elem && _head->_next != NULL)
                {
                    /*
                    if (_head->_elem == NULL || _tail->_elem == NULL || _head->_next == NULL || _tail->_prev == NULL)
                   {
                        std::cout << this << " PB: CORRUPTION ON ONE." << std::endl;
                   }
                    */

                    CASDequeElement<T> * old_last = _tail->_elem;
                    old_last->_next = e;
                    e->_prev = old_last;
                    _tail->_elem = e;
                    _tail->_prev = e->_prev;

//                    std::cout << this << " push-back on more" << std::endl;
//                    print();
                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }

            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }

    } while (! ok);
}

template <typename T>
T CASDeque<T>::pop_front()
{
    T retval(0);
    CASDequeElement<T> * d(NULL);
    bool ok(false);

    do
    {
        d = NULL;
        retval = 0;
        if (_head->_elem == NULL && _tail->_elem == NULL) // assume zero
        {
            ok = true; // or should we check???
        }
        else if (_head->_next == NULL && _tail->_prev == NULL) // assume one
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev == NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == _tail->_elem && _head->_next == NULL)
                {
                    d = _head->_elem;
                    retval = d->_data;

                    _head->_elem = NULL;
                    _tail->_elem = NULL;
//                    std::cout << this << " pop-front on one" << std::endl;
//                    print();
                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }

        }
        else if (_head->_next == _tail->_elem && _tail->_prev == _head->_elem) // assume = 2
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == _tail->_prev && _head->_next == _tail->_elem)
                {
                    d = _head->_elem;
                    retval = d->_data;

                    _head->_elem = _tail->_elem;
                    _head->_elem->_next = NULL;
                    _tail->_elem->_prev = NULL;
                    _head->_next = NULL;
                    _tail->_prev = NULL;
//                    std::cout << this << " pop-front on two" << std::endl;
//                    print();

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
        else if (_head->_next == _tail->_prev) // assume 3
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_prev != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_next == _tail->_prev)
                {
                    d = _head->_elem;
                    retval = d->_data;

                    _head->_elem = _head->_next;
                    _head->_elem->_prev = NULL;
//                    _head->_prev = NULL;
                    _head->_next = _head->_elem->_next;

//                    std::cout << this << " pop-front on three" << std::endl;
//                    print(); 

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;

                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
        else // assume more
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem != NULL && _head->_next != _tail->_elem && _head->_next != _tail->_prev)
                {
                    d = _head->_elem;
                    retval = d->_data;

                    // Set successor to be new head element
                    _head->_elem = d->_next;

                    _head->_elem->_prev = NULL; // it had "d" as predecessor
                    _head->_prev = NULL;        // transfer changes to _head
                    _head->_next = _head->_elem->_next; // transfer changes to _head

//                    std::cout << "pop-front on more" << std::endl;
//                    print();

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
    } while (! ok);

    /*
    if (d != 0)
        std::cout << this << " front: returning " << d << std::endl;
    */
    delete d;
    return retval;
}

template <typename T>
T CASDeque<T>::pop_back()
{
    T retval(0);
    CASDequeElement<T> * d(NULL);
    bool ok(false);

    do
    {
        d = NULL;
        retval = 0;
        if (_head->_elem == NULL && _tail->_elem == NULL) // assume zero
        {
            ok = true; // or should we check???
        }
        else if (_head->_next == NULL && _tail->_prev == NULL) // assume one
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev == NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == _tail->_elem && _head->_next == NULL)
                {
                    d = _head->_elem;
                    retval = d->_data;

                    _head->_elem = NULL;
                    _tail->_elem = NULL;
                    //std::cout << this << " pop-back on one" << std::endl;
                    //print();
                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }

        }
        else if (_head->_next == _tail->_elem && _tail->_prev == _head->_elem) // assume = 2
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem == _tail->_prev && _head->_next == _tail->_elem)
                {
                    d = _tail->_elem;
                    retval = d->_data;
                    _tail->_elem = _head->_elem;
                    _tail->_elem->_prev = NULL;
                    _head->_elem->_next = NULL;
                    _head->_next = NULL;
                    _tail->_prev = NULL;
                    //std::cout << this << " pop-back on two" << std::endl;
                    //print();

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;

                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
        else if (_head->_next == _tail->_prev) // assume 3
        {
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_prev != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_next == _tail->_prev)
                {
                    d = _tail->_elem;
                    retval = d->_data;

                    _tail->_elem = _tail->_prev;
                    _tail->_elem->_next = NULL;
//                    _tail->_next = NULL;
                    _tail->_prev = _tail->_elem->_prev;

                    //std::cout << this << " pop-back on three" << std::endl;
                    //print();

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;

                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }

        }
        else // assume more
        {
            // it must be somehow made sure that there is no NULL/inconsistent case here
            bool blockT = CAS(&_tail->_blocked, 0, 1);
            if (blockT && _tail->_elem != NULL && _tail->_prev != NULL)
            {
                bool blockH = CAS(&_head->_blocked, 0, 1);
                if (blockH && _head->_elem != NULL && _head->_next != _tail->_elem && _head->_next != _tail->_prev)
                {
                    d = _tail->_elem;
                    retval = d->_data;

                        // Set successor to be new head element
                    _tail->_elem = d->_prev;

                    _tail->_elem->_next = NULL; // it had "d" as successor
                    _tail->_next = NULL;        // transfer changes to _head
                    _tail->_prev = _tail->_elem->_prev; // transfer changes to _head

//                    std::cout << "pop-back on more" << std::endl;
//                    print();

                    _head->_blocked = 0;
                    _tail->_blocked = 0;
                    ok = true;
                }
                else
                {
                    if (blockH)
                        _head->_blocked = 0;

                    _tail->_blocked = 0;
                }
            }
            else
            {
                if (blockT)
                    _tail->_blocked = 0;
            }
        }
    } while (! ok);

    /*
    if (d != 0)
        std::cout << this << " back: returning " << d << std::endl;
    */

    delete d;
    return retval;
}

template <typename T>
void CASDeque<T>::print()
{
    std::cout << "List " << this << std::endl;
    std::cout << "head = " << _head << " - elem = " << _head->_elem << " - next = " << _head->_next << std::endl;
    std::cout << "tail = " << _tail << " - elem = " << _tail->_elem << " - prev = " << _tail->_prev << std::endl;

    CASDequeEndElement<T> * h = _head;

    if (h->_elem != NULL)
        std::cout << h->_elem << " -> ";
    else
        return;

    CASDequeElement<T> * cur = _head->_elem;
    while (cur->_next != NULL)
    {
        cur = cur->_next;
        std::cout << cur << " -> ";
    }
    std::cout << std::endl;

}
