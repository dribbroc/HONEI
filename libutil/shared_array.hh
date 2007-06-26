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

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/log.hh>
#include <libutil/stringify.hh>

#include <map>
#include <iostream>

namespace pg512 ///< \todo Namespace name?
{
    class SharedArrayError :
        public Exception
    {
        public:
            SharedArrayError(const std::string & message) throw();
    };

    template <typename DataType_> class SharedArrayCounter
    {
        private:
            SharedArrayCounter()
            {
            }

            /// Our map of reference counters.
            std::map<DataType_ *, signed long> _counter_map;

        public:
            static SharedArrayCounter * instance()
            {
                static SharedArrayCounter result;

                return &result;
            }

            void increase(DataType_ * entry)
            {
                if (_counter_map.find(entry) == _counter_map.end())
                    _counter_map[entry] = 0;

                ++_counter_map[entry];
            }

            void decrease(DataType_ * entry)
            {
                if (_counter_map.find(entry) == _counter_map.end())
                    _counter_map[entry] = 0;

                --_counter_map[entry];
            }

            signed long count(DataType_ * entry)
            {
                if (_counter_map.find(entry) == _counter_map.end())
                    return 0;
                else
                    return _counter_map[entry];
            }
    };

    template <typename DataType_> class SharedArray
    {
        private:
            /// Pointer to our array.
            mutable DataType_ *_array;

        public:
            /// Constructor.
            SharedArray(DataType_ *array) :
                _array(array)
            {
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) >= 0, "reference counter for array '" +
                    stringify(array) + "' is below zero!");
                SharedArrayCounter<DataType_>::instance()->increase(_array);
            }

            /// Copy-Constructor.
            SharedArray(const SharedArray<DataType_> & other) :
                _array(other._array)
            {
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) >= 0, "reference counter for array '" +
                    stringify(other._array) + "' is below zero!");
                SharedArrayCounter<DataType_>::instance()->increase(_array);
            }

            /// Destructor.
            ~SharedArray()
            {
                SharedArrayCounter<DataType_>::instance()->decrease(_array);

                signed long count(SharedArrayCounter<DataType_>::instance()->count(_array));
                ASSERT(count >= 0, "reference counter for array '" + stringify(_array) + "' is below zero!");

                if (count == 0)
                    delete[] _array;
            }

            DataType_ & operator[] (std::ptrdiff_t index) const
            {
                ASSERT(index >= 0, "index '" + stringify(index) + "' is out of bounds!");
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

            void reset(DataType_ * array) const
            {
                if (_array != array)
                    Log::instance()->message(ll_full, "reset: new array = " + stringify(array) + " taken with count = "
                            + stringify(SharedArrayCounter<DataType_>::instance()->count(array)));

                // Increase first to avoid race conditions like 'A.reset(A.get())'
                SharedArrayCounter<DataType_>::instance()->increase(array);
                SharedArrayCounter<DataType_>::instance()->decrease(_array);

                signed long count(SharedArrayCounter<DataType_>::instance()->count(_array));
                ASSERT(count >= 0, "reference counter for array '" + stringify(_array) + "' is below zero!");

                if (count == 0)
                    delete[] _array;

                _array = array;
            }
    };
}

#endif
