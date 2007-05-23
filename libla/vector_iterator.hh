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

#include <libla/vector.hh>

#include <iterator>
#include <tr1/memory>

namespace pg512 ///< \todo Namespace name?
{
    template <typename DataType_, typename VectorType_ = Vector<DataType_> > class VectorIteratorBase :
        public std::iterator<std::forward_iterator_tag, VectorType_>
    {
        public:
            /// Preincrement operator.
            virtual VectorIteratorBase<DataType_, VectorType_> & operator++ () = 0;

            /// Returns our index.
            virtual const unsigned long index() const = 0;
    };

    template <typename DataType_, typename VectorType_ = Vector<DataType_> > class VectorIteratorImplBase :
        public VectorIteratorBase<DataType_, VectorType_>
    {
        public:
            /// Equality operator.
            virtual bool operator== (const VectorIteratorImplBase<DataType_, VectorType_> & other) const = 0;

            /// Inequality operator.
            virtual bool operator!= (const VectorIteratorImplBase<DataType_, VectorType_> & other) const = 0;

            /// Dereference operator.
            virtual VectorType_ & operator* () const = 0;

            /// Returns our parent.
            virtual const Matrix<DataType_> * parent() const = 0;
    };

    template <typename DataType_, typename VectorType_ = Vector<DataType_> > class VectorIteratorWrapper :
        public VectorIteratorBase<DataType_, VectorType_>
    {
        private:
            /// Our wrapped iterator.
            std::tr1::shared_ptr<VectorIteratorImplBase<DataType_, VectorType_> > _iterator;

        public:
            /// Constructor.
            VectorIteratorWrapper(VectorIteratorImplBase<DataType_, VectorType_> *iterator) :
                _iterator(iterator)
            {
            }

            /// Copy-constructor.
            VectorIteratorWrapper(const VectorIteratorWrapper<DataType_, VectorType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Preincrement operator.
            virtual VectorIteratorBase<DataType_, VectorType_> & operator++ ()
            {
                ++(*_iterator);
                return *this;
            }

            /// Postincrement operator.
            virtual VectorIteratorWrapper<DataType_, VectorType_> operator++ (int)
            {
                VectorIteratorWrapper<DataType_, VectorType_> result(*this);

                ++(*_iterator);

                return result;
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

            /// Dereference operator 
            virtual VectorType_ & operator* () const
            {
                return **_iterator;
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _iterator->index();
            }
    };
}

#endif
