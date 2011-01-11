/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#pragma once
#ifndef LIBUTIL_GUARD_SHARED_ARRAY_IMPL_HH
#define LIBUTIL_GUARD_SHARED_ARRAY_IMPL_HH 1

#include <honei/util/assertion.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array.hh>
#include <honei/util/stringify.hh>
#include <honei/util/type_traits.hh>

#include <algorithm>

namespace honei
{
    template <typename DataType_> struct Implementation<SharedArray<DataType_> >
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
            Lock l(*mutex);
            TypeTraits<DataType_>::create(array, s, DataType_());
        }

        /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
        Implementation(const Implementation & other);

        /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
        Implementation & operator= (const Implementation & other);

        /// Destructor.
        ~Implementation()
        {
            {
                Lock l(*mutex);
                TypeTraits<DataType_>::destroy(array, size);
                TypeTraits<DataType_>::free(array, size);
            }

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
        PrivateImplementationPattern<SharedArray<DataType_>, Shared>(new Implementation<SharedArray<DataType_> >(size))
    {
    }

    template <typename DataType_>
    SharedArray<DataType_>::SharedArray(const SharedArray<DataType_> & other) :
        PrivateImplementationPattern<SharedArray<DataType_>, Shared>(other._imp)
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
                stringify(this->_imp->size) + "':");
        //ASSERT(index >= 0, "index '" + stringify(index) + "' is out of bounds!");
        ASSERT(index < this->_imp->size, "index '" + stringify(index) + "' is out of bounds!");

        return this->_imp->array[index];
    }

    template <typename DataType_>
    bool SharedArray<DataType_>::operator! () const
    {
        return this->_imp->array == 0;
    }

    template <typename DataType_>
    unsigned long SharedArray<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    DataType_ * SharedArray<DataType_>::get() const
    {
        return this->_imp->array;
    }

    template <typename DataType_>
    void SharedArray<DataType_>::reset(unsigned long size, DataType_ * array) const
    {
        CONTEXT("When resetting SharedArray of size '" + stringify(this->_imp->size) + "' with POA '" +
                stringify(array) + "':");
        ASSERT(this->_imp->array != array, "new array is identical with old array!");
        Lock l(*this->_imp->mutex);

        this->_imp->reset(size, array);
    }
}

#endif
