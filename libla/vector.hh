/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_VECTOR_HH
#define LIBLA_GUARD_VECTOR_HH 1

#include <libla/element_iterator.hh>
#include <libla/vector_error.hh>
#include <libutil/shared_array.hh>

#include <cmath>
#include <iterator>
#include <limits>
#include <ostream>

/**
 * \file
 *
 * Interface declarations for Vector-based types.
 *
 * \ingroup grpvector
 **/
namespace honei
{
    /**
     * \brief Vector is the abstract baseclass for all vector-like types used.
     *
     * \ingroup grpvector
     **/
    template <typename DataType_> class Vector
    {
        protected:
            template <typename ElementType_> class ElementIteratorBase;
            typedef ElementIteratorBase<DataType_> VectorElementIterator;

        public:
            /// Type of the const iterator over our elements.
            template <typename ElementType_> class ConstElementIteratorWrapper;
            typedef ConstElementIteratorWrapper<DataType_> ConstElementIterator;

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const = 0;

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const = 0;

            /// Returns const iterator pointing to a given element of the vector.
            virtual ConstElementIterator element_at(unsigned long index) const = 0;

            /// Type of the iterator over our elements.
            template <typename ElementType_> class ElementIteratorWrapper;
            typedef ElementIteratorWrapper<DataType_> ElementIterator;

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements() = 0;

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() = 0;

            /// Returns const iterator pointing to a given element of the vector.
            virtual ElementIterator element_at(unsigned long index) = 0;

            /// Returns our size.
            virtual unsigned long size() const = 0;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const = 0;

            /// Retrieves element by index, zero-based, assignable
            virtual DataType_ & operator[] (unsigned long index) = 0;

            /// Returns a pointer to a copy.
            virtual Vector * copy() const = 0;
    };

    /**
     * \brief Vector::ElementIteratorBase declares the minimal interface for any ElementIterator implementation
     * \brief for vector-like types.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class Vector<DataType_>::ElementIteratorBase :
        public IteratorBase<DataType_, Vector<DataType_> >
    {
    };

    /**
     * \brief Vector::ElementIteratorWrapper provides a covariant mutable forward iterator that wraps the actual
     * \brief ElementIteratorBase implementations of any of Vector's descendants.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class Vector<DataType_>::ElementIteratorWrapper<DataType_> :
        public std::iterator<std::forward_iterator_tag, DataType_>
    {
        private:
            std::tr1::shared_ptr<ElementIteratorBase<DataType_> > _iterator;

        public:
            friend class ConstElementIteratorWrapper<DataType_>;

            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param iterator An instance of one of ElementIteratorBase's descendants that shall be wrapped.
             */
            ElementIteratorWrapper(ElementIteratorBase<DataType_> * iterator) :
                _iterator(iterator)
            {
                if (! iterator)
                    throw std::string("Eek. Iterator is 0, that should not happen!");
            }

            /// Copy-constructor.
            ElementIteratorWrapper(const ElementIteratorWrapper<DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual ElementIteratorWrapper<DataType_> & operator++ ()
            {
                ++(*_iterator);

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                return *(*_iterator);
            }

            /// Comparison operator for equality.
            virtual bool operator== (const ElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const ElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns our index.
            unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns a pointer to our parent container.
            const Vector<DataType_> * parent() const
            {
                return _iterator->parent();
            }

            /// \}
    };

    /**
     * \brief Vector::ConstElementIteratorWrapper provides a covariant const forward iterator that wraps the actual
     * \brief ElementIteratorBase implementations of any of Vector's descendants.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class Vector<DataType_>::ConstElementIteratorWrapper<DataType_> :
        public std::iterator<std::forward_iterator_tag, const DataType_>
    {
        private:
            std::tr1::shared_ptr<ElementIteratorBase<DataType_> > _iterator;

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param iterator An instance of one of ElementIteratorBase's descendants that shall be wrapped.
             */
            ConstElementIteratorWrapper(ElementIteratorBase<DataType_> * iterator) :
                _iterator(iterator)
            {
                if (! iterator)
                    throw std::string("Eek. Iterator is 0, that should not happen!");
            }

            /// Copy-constructor.
            ConstElementIteratorWrapper(const ConstElementIteratorWrapper<DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Conversion-constructor from mutable ElementIteratorWrapper.
            ConstElementIteratorWrapper(const ElementIteratorWrapper<DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual ConstElementIteratorWrapper<DataType_> & operator++ ()
            {
                ++(*_iterator);

                return *this;
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                const ElementIteratorBase<DataType_> & iterator(*_iterator);

                return *iterator;
            }

            /// Comparison operator for equality.
            bool operator== (const ConstElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Comparison operator for inequality.
            bool operator!= (const ConstElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns our index.
            unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns a pointer to our parent container.
            const Vector<DataType_> * parent() const
            {
                return _iterator->parent();
            }

            /// \}
    };

    /**
     * \brief Compare two Vectors for equality.
     *
     * \ingroup grpvectoroperations
     */
    template <typename DataType_> bool operator== (const Vector<DataType_> & left, const Vector<DataType_> & right)
    {
        bool result(true);
        if (left.size() != right.size())
            throw VectorSizeDoesNotMatch(left.size(), right.size());

        for (typename Vector<DataType_>::ConstElementIterator i(left.begin_elements()), i_end(left.end_elements()),
                j(right.begin_elements()) ; i != i_end ; ++i)
        {
            if (fabs(*i - *j) <= std::numeric_limits<DataType_>::epsilon())
            {
                ++j;
                continue;
            }

            result = false;
            break;
        }

        return result;
    }

    /**
     * \brief Output our Vector to an ostream.
     *
     * \ingroup grpvectoroperations
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const Vector<DataType_> & v)
    {
        lhs << "[ ";
        for (typename Vector<DataType_>::ConstElementIterator i(v.begin_elements()), i_end(v.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << *i << " ";
        }
        lhs << "]";

        return lhs;
    }
}

#endif
