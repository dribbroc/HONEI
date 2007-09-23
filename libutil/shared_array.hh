/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_SHARED_ARRAY_HH
#define LIBUTIL_GUARD_SHARED_ARRAY_HH 1

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/log.hh>
#include <libutil/stringify.hh>

#include <map>
#include <iostream>

namespace honei
{
    /**
     * SharedArrayError is thrown by SharedArray and related classes.
     *
     * \ingroup grpexceptions
     * \ingroup grpsharedarray
     */
    class SharedArrayError :
        public Exception
    {
        public:
            SharedArrayError(const std::string & message) throw();
    };

    /**
     * SharedArrayCounter is the counter class template for SharedArray.
     *
     * \warning For internal use by SharedArray only!
     */
    template <typename DataType_> class SharedArrayCounter
    {
        private:
            /// \name Constructors
            /// \{

            /// Constructor.
            SharedArrayCounter()
            {
            }

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            SharedArrayCounter(const SharedArrayCounter<DataType_> &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            SharedArrayCounter<DataType_> & operator= (const SharedArrayCounter<DataType_> &);

            /// \}

            /// Our map of reference counters.
            std::map<DataType_ *, signed long> _counter_map;

        public:
            /// Return the only instance of SharedArrayCounter.
            static SharedArrayCounter<DataType_> * instance()
            {
                static SharedArrayCounter<DataType_> result;

                return &result;
            }

            /// Increase counter for given pointer.
            void increase(DataType_ * entry)
            {
                CONTEXT("When increasing counter for entry '" + stringify(entry) + "':");

                if (_counter_map.find(entry) == _counter_map.end())
                    _counter_map[entry] = 0;

                ++_counter_map[entry];
            }

            /// Decrease counter for given pointer.
            void decrease(DataType_ * entry)
            {
                CONTEXT("When decreasing counter for entry '" + stringify(entry) + "':");

                if (_counter_map.find(entry) == _counter_map.end())
                    _counter_map[entry] = 0;

                --_counter_map[entry];
            }

            signed long count(DataType_ * entry)
            {
                CONTEXT("When retrieving counter value for entry '" + stringify(entry) + "':");

                if (_counter_map.find(entry) == _counter_map.end())
                    return 0;
                else
                    return _counter_map[entry];
            }
    };

    /**
     * SharedArray is a class template that allows sharing one array among several owners.
     */
    template <typename DataType_> class SharedArray
    {
        private:
            /// Our data.
            mutable DataType_ *_array;

            /// Our size.
            mutable unsigned long _size;

        public:
            /// \name Constructors and destructor.
            /// \{

            /// (Explicit) constructor.
            explicit SharedArray(unsigned long size) :
                _array(new DataType_[size]),
                _size(size)
            {
                CONTEXT("When creating SharedArray of size '" + stringify(size) + "' with default elements:");
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) == 0, "reference counter for array '" +
                    stringify(_array) + "' is not zero!");

                SharedArrayCounter<DataType_>::instance()->increase(_array);
            }

            /// Constructor.
            SharedArray(unsigned long size, DataType_ * array) :
                _array(array),
                _size(size)
            {
                CONTEXT("When creating SharedArray of size '" + stringify(size) + "' with given POA:");
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) >= 0, "reference counter for array '" +
                    stringify(_array) + "' is below zero!");

                SharedArrayCounter<DataType_>::instance()->increase(_array);
            }

            /// Copy-constructor.
            SharedArray(const SharedArray<DataType_> & other) :
                _array(other._array),
                _size(other._size)
            {
                CONTEXT("When creating SharedArray of size '" + stringify(_size) + "' from given SharedArray:");
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) >= 0, "reference counter for array '" +
                    stringify(other._array) + "' is below zero!");

                SharedArrayCounter<DataType_>::instance()->increase(_array);
            }

            /// Assignment-operator.
            SharedArray<DataType_> & operator= (const SharedArray<DataType_> & other)
            {
                CONTEXT("When assigning SharedArray to SharedArray:");
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(_array) >= 0, "reference counter for other's "
                        "array is below zero!");

                // Increase first to avoid race conditions like 'A.reset(A.get())'
                SharedArrayCounter<DataType_>::instance()->increase(other._array);
                SharedArrayCounter<DataType_>::instance()->decrease(_array);

                _array = other._array;
                _size = other._size;
            }

            /// Destructor.
            ~SharedArray()
            {
                CONTEXT("When destroying SharedArray of size '" + stringify(_size) + "':");

                SharedArrayCounter<DataType_>::instance()->decrease(_array);

                signed long count(SharedArrayCounter<DataType_>::instance()->count(_array));
                ASSERT(count >= 0, "reference counter for array '" + stringify(_array) + "' is below zero!");

                if (count == 0)
                    delete[] _array;
            }

            /// \}

            /// \name Operators to mimic POA
            /// \{

            /// Subscript operator, return element at given index.
            DataType_ & operator[] (std::ptrdiff_t index) const
            {
                CONTEXT("When accessing SharedArray-element at index '" + stringify(index) + "':");
                ASSERT(index >= 0, "index '" + stringify(index) + "' is out of bounds!");
                ASSERT(index < _size, "index '" + stringify(index) + "' is out of bounds!");

                return _array[index];
            }

            /// Return whether we hold any data.
            inline bool operator! () const
            {
                return _array == 0;
            }

            /// \}

            /// Return our size.
            inline unsigned long size() const
            {
                return _size;
            }

            /// Return our POA.
            inline DataType_ * get() const
            {
                return _array;
            }

            /// Reset us with a new size and a new POA.
            void reset(unsigned long size, DataType_ * array) const
            {
                CONTEXT("When resetting SharedArray of size '" + stringify(_size) + "' with POA '" +
                        stringify(array) + "':");
                ASSERT(_array != array, "new array is identical with old array!");
                ASSERT(SharedArrayCounter<DataType_>::instance()->count(array) == 0,
                        "reference counter for new array '" + stringify(array) + "' is not zero!");

                // Increase first to avoid race conditions like 'A.reset(A.get())'
                SharedArrayCounter<DataType_>::instance()->increase(array);
                SharedArrayCounter<DataType_>::instance()->decrease(_array);

                signed long count(SharedArrayCounter<DataType_>::instance()->count(_array));
                ASSERT(count >= 0, "reference counter for new array '" + stringify(_array) + "' is below zero!");

                if (count == 0)
                    delete[] _array;

                _array = array;
                _size = size;
            }
    };
}

#endif
