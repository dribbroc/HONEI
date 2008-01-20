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

#ifndef LIBLA_GUARD_VECTOR_ITERATOR_HH
#define LIBLA_GUARD_VECTOR_ITERATOR_HH 1

#include <honei/libla/matrix.hh>

#include <iterator>
#include <tr1/memory>

namespace honei
{
    class VectorIteratorTraits
    {
        public:
            /// Returns true if the referenced vector already exists.
            virtual bool exists() const = 0;

            /// Returns our index.
            virtual const unsigned long index() const = 0;

            /// Destructor.
            virtual ~VectorIteratorTraits() {}
    };

    template <typename DataType_, typename VectorType_> class VectorIteratorBase :
        public std::iterator<std::forward_iterator_tag, VectorType_>,
        public VectorIteratorTraits
    {
        public:
            /// Preincrement operator.
            virtual VectorIteratorBase<DataType_, VectorType_> & operator++ () = 0;

            /// Equality operator.
            virtual bool operator== (const VectorIteratorBase<DataType_, VectorType_> &) const = 0;

            /// Inequality operator.
            virtual bool operator!= (const VectorIteratorBase<DataType_, VectorType_> &) const = 0;

            /// Dereference operator (mutable).
            virtual VectorType_ & operator* () = 0;

            /// Dereference operator (const).
            virtual const VectorType_ & operator* () const = 0;

            /// Returns our parent.
            virtual const Matrix<DataType_> * parent() const = 0;
    };

    template <typename DataType_, typename VectorType_> class ConstVectorIteratorWrapper;

    template <typename DataType_, typename VectorType_> class VectorIteratorWrapper :
        public std::iterator<std::forward_iterator_tag, VectorType_>,
        public VectorIteratorTraits
    {
        private:
            /// Our wrapped iterator.
            std::tr1::shared_ptr<VectorIteratorBase<DataType_, VectorType_> > _iterator;

        public:
            friend class ConstVectorIteratorWrapper<DataType_, VectorType_>;

            /// Constructor.
            VectorIteratorWrapper(VectorIteratorBase<DataType_, VectorType_> *iterator) :
                _iterator(iterator)
            {
            }

            /// Copy-constructor.
            VectorIteratorWrapper(const VectorIteratorWrapper<DataType_, VectorType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Preincrement operator.
            virtual VectorIteratorWrapper<DataType_, VectorType_> & operator++ ()
            {
                ++(*_iterator);
                return *this;
            }

            /// Equality operator.
            virtual bool operator== (const VectorIteratorWrapper<DataType_, VectorType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Inequality operator.
            virtual bool operator!= (const VectorIteratorWrapper<DataType_, VectorType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// Dereference operator (return reference to type).
            virtual VectorType_ & operator* ()
            {
                return **_iterator;
            }

            /// Dereference operator (return pointer to type).
            virtual VectorType_ * operator-> ()
            {
                return &(*(*_iterator));
            }

            /// Returns true if the referenced vector already exists.
            virtual bool exists() const
            {
                return _iterator->exists();
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns a pointer to our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return _iterator->parent();
            }
    };

    template <typename DataType_, typename VectorType_> class ConstVectorIteratorWrapper :
        public std::iterator<std::forward_iterator_tag, VectorType_>,
        public VectorIteratorTraits
    {
        private:
            /// Our wrapped iterator.
            std::tr1::shared_ptr<VectorIteratorBase<DataType_, VectorType_> > _iterator;

        public:
            /// Constructor.
            ConstVectorIteratorWrapper(VectorIteratorBase<DataType_, VectorType_> *iterator) :
                _iterator(iterator)
            {
            }

            /// Conversion-constructor.
            ConstVectorIteratorWrapper(const VectorIteratorWrapper<DataType_, VectorType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Copy-constructor.
            ConstVectorIteratorWrapper(const ConstVectorIteratorWrapper<DataType_, VectorType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Preincrement operator.
            virtual ConstVectorIteratorWrapper<DataType_, VectorType_> & operator++ ()
            {
                ++(*_iterator);
                return *this;
            }

            /// Equality operator.
            virtual bool operator== (const ConstVectorIteratorWrapper<DataType_, VectorType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Inequality operator.
            virtual bool operator!= (const ConstVectorIteratorWrapper<DataType_, VectorType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// Dereference operator (return reference to type).
            virtual const VectorType_ & operator* () const
            {
                const VectorIteratorBase<DataType_, VectorType_> & iterator(*_iterator);

                return *iterator;
            }

            /// Dereference operator (return pointer to type).
            virtual const VectorType_ * operator-> () const
            {
                const VectorIteratorBase<DataType_, VectorType_> & iterator(*_iterator);

                return &(*iterator);
            }

            /// Returns true if the referenced vector already exists.
            virtual bool exists() const
            {
                return _iterator->exists();
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns a pointer to our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return _iterator->parent();
            }
    };
}

#endif
