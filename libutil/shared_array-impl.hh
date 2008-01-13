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

#ifndef LIBUTIL_GUARD_SHARED_ARRAY_IMPL_HH
#define LIBUTIL_GUARD_SHARED_ARRAY_IMPL_HH 1

#include <libutil/shared_array.hh>

#include <libutil/assertion.hh>
#include <libutil/lock.hh>
#include <libutil/mutex.hh>
#include <libutil/stringify.hh>
#include <libutil/type_traits.hh>

#include <algorithm>

namespace honei
{
    template <typename DataType_> struct SharedArray<DataType_>::Implementation
    {
        /// Our data.
        DataType_ * array;

        /// Our size.
        unsigned long size;

        /// Our mutex.
        Mutex * const mutex;

        /// \name Basic operations
        /// \{

        /// Constructor.
        Implementation(unsigned long s) :
            array(TypeTraits<DataType_>::allocate(s)), // Needs to be 'free'ed!
            size(s),
            mutex(new Mutex)
        {
            TypeTraits<DataType_>::create(array, s, DataType_());
        }

        /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
        Implementation(const Implementation & other);

        /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
        Implementation & operator= (const Implementation & other);

        /// Destructor.
        ~Implementation()
        {
            TypeTraits<DataType_>::destroy(array, size);
            TypeTraits<DataType_>::free(array, size);

            delete mutex;
        }

        /// \}

        /// Reset us with a new array.
        inline void reset(unsigned long s, DataType_ * a)
        {
            TypeTraits<DataType_>::destroy(array, size);
            TypeTraits<DataType_>::free(array, size);

            array = a;
            size = s;
        }
    };

    template <typename DataType_>
    SharedArray<DataType_>::SharedArray(unsigned long size) :
        _imp(new Implementation(size))
    {
    }

    template <typename DataType_>
    SharedArray<DataType_>::SharedArray(const SharedArray<DataType_> & other) :
        _imp(other._imp)
    {
    }

    template <typename DataType_>
    SharedArray<DataType_>::~SharedArray()
    {
    }

    template <typename DataType_>
    DataType_ & SharedArray<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When accessing SharedArray-element at index '" + stringify(index) + "' in array of size '" +
                stringify(_imp->size) + "':");
        ASSERT(index >= 0, "index '" + stringify(index) + "' is out of bounds!");
        ASSERT(index < _imp->size, "index '" + stringify(index) + "' is out of bounds!");

        return _imp->array[index];
    }

    template <typename DataType_>
    bool SharedArray<DataType_>::operator! () const
    {
        return _imp->array == 0;
    }

    template <typename DataType_>
    unsigned long SharedArray<DataType_>::size() const
    {
        return _imp->size;
    }

    template <typename DataType_>
    DataType_ * SharedArray<DataType_>::get() const
    {
        return _imp->array;
    }

    template <typename DataType_>
    void SharedArray<DataType_>::reset(unsigned long size, DataType_ * array) const
    {
        CONTEXT("When resetting SharedArray of size '" + stringify(_imp->size) + "' with POA '" +
                stringify(array) + "':");
        ASSERT(_imp->array != array, "new array is identical with old array!");
        Lock l(*_imp->mutex);

        _imp->reset(size, array);
    }
}

#endif
