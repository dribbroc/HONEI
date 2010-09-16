/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/segment_list.hh>
#include <honei/util/lock.hh>

#include <tr1/memory>

using namespace honei;
using namespace honei::mc;

//#define LIST_DEBUG2
//#define SIZE_DEBUG

SegmentList::SegmentList(const unsigned procs) :
    _num_segs(procs),
    _size(0),
    _segs(new Segment*[procs]),
    _global_mutex(new Mutex)
{
    for (unsigned i(0) ; i < procs ; ++i)
        _segs[i] = new Segment;
}

SegmentList::~SegmentList()
{
    for (unsigned i(0) ; i < _num_segs ; ++i)
        delete _segs[i];

    delete _global_mutex;

    delete[] _segs;
}

void SegmentList::push_back(ThreadTask * elem)
{
    unsigned seg_id(0);

    while (_segs[seg_id]->count > (1 + (_size / _num_segs)))
        ++seg_id;

    _segs[seg_id]->lock();

#ifdef LIST_DEBUG
    std::cout << "Pushing in Segment : " << seg_id << " with size = " << _segs[seg_id]->count << std::endl;
#endif

    _segs[seg_id]->push(elem);

    {
        Lock l(*_global_mutex);
        ++_size;
    }

#ifdef LIST_DEBUG2
    std::cout << "Pushed " << elem << " in Segment : " << seg_id << " elem with sched_id " << elem->ticket->sid_min() << std::endl;
//    std::cout << "Will now unlock. size = " << _segs[seg_id]->count << std::endl;
#endif

    _segs[seg_id]->unlock();

#ifdef SIZE_DEBUG
    {
//        Lock l(*_global_mutex);
        std::cout << "Leaving push with size = " << _size << std::endl;
    }
#endif

}

ThreadTask * SegmentList::extract(std::tr1::function<bool (ThreadTask * const)> & comp)
{
    ThreadTask * res(NULL);
    unsigned seg_id(0);

    while (seg_id != _num_segs)
    {
        _segs[seg_id]->lock();

#ifdef LIST_DEBUG
    std::cout << "Thread " << (comp.target<TaskComp>())->sched_lpu << " Searching in Segment : " << seg_id << std::endl;
#endif

        res = _segs[seg_id]->find_and_extract(comp);

        if (res != NULL)
        {
            {
                Lock l(*_global_mutex);
                --_size;
            }
#ifdef LIST_DEBUG
    std::cout << "Found elem in Segment : " << seg_id << std::endl;
    std::cout << "Will now unlock. size = " << _segs[seg_id]->count << std::endl;
#endif
            _segs[seg_id]->unlock();

            break;
        }
        else
        {

#ifdef LIST_DEBUG2
    std::cout << "Thread " << (comp.target<TaskComp>())->sched_lpu << " found nothing in segment : " << seg_id << " with size " << _segs[seg_id]->count << std::endl;
#endif

#ifdef LIST_DEBUG
    std::cout << "Found nothing in Segment : " << seg_id << std::endl;
    std::cout << "Will now unlock. size = " << _segs[seg_id]->count << std::endl;
    if (_segs[seg_id]->head != NULL) std::cout << " head elem = " << _segs[seg_id]->head->data << std::endl;
#endif
            _segs[seg_id]->unlock();

            ++seg_id;
        }
    }

#ifdef SIZE_DEBUG
    {
//        Lock l(*_global_mutex);
        std::cout << "Leaving extract with size = " << _size << std::endl;
    }
#endif

    return res;
}
