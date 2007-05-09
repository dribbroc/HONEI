/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBUTIL_GUARD_SHARED_ARRAY_HH
#define LIBUTIL_GUARD_SHARED_ARRAY_HH 1

#include "exception.hh"

#include <map>

namespace pg512 ///< \todo Namespace name?
{
    class SharedArrayError :
        public Exception
    {
        public:
            SharedArrayError(const std::string & message) throw();
    };

    template <typename DataType_> class SharedArray
    {
        private:
            /// Pointer to our array.
            DataType_ *_array;

            /// Our map of reference counters.
            static std::map<DataType_ *, unsigned long> _counter_map;

        public:
            /// Constructor.
            SharedArray(DataType_ *array) :
                _array(array)
            {
                _counter_map[_array]++;
            }

            /// Constructor.
            SharedArray(const SharedArray<DataType_> & other) :
                _array(other._array)
            {
                _counter_map[_array]++;
            }

            /// Destructor.
            ~SharedArray()
            {
                _counter_map[_array]--;
#if 0
                if (_counter_map[_array] < 0)
                    Log(ll_critical, "~SharedArray: reference counter below zero for address " + stringify(_array));
                else
#endif
                if (_counter_map[_array] == 0)
                    delete[] _array;
            }

            DataType_ & operator* ()
            {
                if (_counter_map[_array] <= 0)
                    throw std::string("ShareArray::operator*(): reference counter below zero for address " + stringify(_array));
                else
                    return _array[0];
            }

            DataType_ & operator[] (std::ptrdiff_t index) const
            {
                if (_counter_map[_array] <= 0)
                    throw std::string("ShareArray::operator*(): reference counter below zero for address " + stringify(_array));
                else
                    return _array[index];
            }

            bool operator! () const
            {
                return _array == 0;
            }

            DataType_ * get() const
            {
                return _array;
            }
    };
}

#endif
