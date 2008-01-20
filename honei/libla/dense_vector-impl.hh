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

#ifndef LIBLA_GUARD_DENSE_VECTOR_IMPL_HH
#define LIBLA_GUARD_DENSE_VECTOR_IMPL_HH 1

#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_vector_range.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/element_iterator.hh>
#include <honei/libutil/assertion.hh>
#include <honei/libutil/shared_array.hh>
#include <honei/libutil/stringify.hh>

#include <algorithm>
#include <string>

namespace honei
{
    // Forward declarations.
    template <typename DataType_> class DenseMatrix;
    template <typename DataType_> class DenseVectorRange;

    /**
     * \brief DenseVector::Implementation is the private implementation class for DenseVector.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVector<DataType_>::Implementation
    {
        private:
            /// Unwanted copy-constructor: Do not implement. See EffC++, Item 27.
            Implementation(const Implementation &);

            /// Unwanted assignment operator: Do not implement. See EffC++, Item 27.
            Implementation & operator= (const Implementation &);

        public:
            /// Our elements.
            SharedArray<DataType_> elements;

            /// Our size.
            const unsigned long size;

            /// Our offset.
            const unsigned long offset;

            /// Our stepsize.
            const unsigned long stepsize;

            /// Constructor.
            Implementation(unsigned long s, unsigned long o, unsigned long ss) :
                elements(o + ss * (s + 1) - 1),
                size(s),
                offset(o),
                stepsize(ss)
            {
            }

            /// Constructor.
            Implementation(const SharedArray<DataType_> & e, unsigned long s, unsigned long o, unsigned long ss) :
                elements(e),
                size(s),
                offset(o),
                stepsize(ss)
            {
            }
    };

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, const SharedArray<DataType_> & elements,
            unsigned long offset, unsigned stepsize) :
        _imp(new Implementation(elements, size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, const unsigned long offset,
            const unsigned long stepsize) :
        _imp(new Implementation(size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, DataType_ value, unsigned long offset,
            unsigned long stepsize) :
        _imp(new Implementation(size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");

        /// \todo Replace by the following once DenseVector has lost Range and Slice functionality:
        // TypeTraits<DataType_>::fill(_imp->elements.get(), size, value);

        DataType_ *  target(_imp->elements.get());
        for (unsigned long i(offset) ; i < (stepsize * size + offset) ; i += stepsize)
            target[i] = value;
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const SparseVector<DataType_> & other) :
        _imp(new Implementation(other.size(), 0, 1))
    {
        CONTEXT("When creating DenseVector form SparseVector:");
        ASSERT(other.size() > 0, "size is zero!");

        TypeTraits<DataType_>::fill(_imp->elements.get(), other.size(), DataType_(0));

        for (typename Vector<DataType_>::ConstElementIterator i(other.begin_non_zero_elements()),
                i_end(other.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            (*this)[i.index()] = *i;
        }
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const DenseVector<DataType_> & other) :
        _imp(new Implementation(other._imp->elements, other._imp->size, other._imp->offset, other._imp->stepsize))
    {
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::begin_elements() const
    {
        return ConstElementIterator(new DenseElementIterator<DataType_>(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::end_elements() const
    {
        return ConstElementIterator(new DenseElementIterator<DataType_>(*this, _imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::element_at(unsigned long index) const
    {
        return ConstElementIterator(new DenseElementIterator<DataType_>(*this, index));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::begin_elements()
    {
        return ElementIterator(new DenseElementIterator<DataType_>(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::end_elements()
    {
        return ElementIterator(new DenseElementIterator<DataType_>(*this, _imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::element_at(unsigned long index)
    {
        return ElementIterator(new DenseElementIterator<DataType_>(*this, index));
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::size() const
    {
        return _imp->size;
    }

    template <typename DataType_>
    const DataType_ & DenseVector<DataType_>::operator[] (unsigned long index) const
    {
        return _imp->elements[_imp->stepsize * index + _imp->offset];
    }

    template <typename DataType_>
    DataType_ & DenseVector<DataType_>::operator[] (unsigned long index)
    {
        return _imp->elements[_imp->stepsize * index + _imp->offset];
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::offset() const
    {
        return _imp->offset;
    }

    template <typename DataType_>
    DenseVectorRange<DataType_> DenseVector<DataType_>::range(unsigned long size, unsigned long offset) const
    {
        DenseVectorRange<DataType_> result((*this) , size, offset);
        return result;
    }

    template <typename DataType_>
    inline DataType_ * DenseVector<DataType_>::elements() const
    {
        return _imp->elements.get();
    }

    template <typename DataType_>
    DenseVector<DataType_> DenseVector<DataType_>::copy() const
    {
        DenseVector result(_imp->size);

        TypeTraits<DataType_>::copy( _imp->elements.get(), result._imp->elements.get(), _imp->size);

        return result;
    }

    /**
     * \brief DenseVector::DenseElementIterator is a simple iterator implementation for dense vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::DenseElementIterator<DataType_> :
        public VectorElementIterator
        {
            private:
                /// Our parent vector.
                const DenseVector<DataType_> & _vector;

                /// Our index.
                unsigned long _index;

            public:
                /// \name Constructors
                /// \{

                /**
                 * Constructor.
                 *
                 * \param vector The parent vector that is referenced by the iterator.
                 * \param index The index into the vector.
                 **/
                DenseElementIterator(const DenseVector<DataType_> & vector, unsigned long index) :
                    _vector(vector),
                    _index(index)
                {
                }

                /// Copy-constructor.
                DenseElementIterator(DenseElementIterator<DataType_> const & other) :
                    _vector(other._vector),
                    _index(other._index)
                {
                }

                /// Destructor.
                virtual ~DenseElementIterator()
                {
                }

                /// \}

                /// \name Forward iterator interface
                /// \{

                /// Preincrement operator.
                virtual DenseElementIterator<DataType_> & operator++ ()
                {
                    CONTEXT("When incrementing iterator by one:");

                    ++_index;

                    return *this;
                }

                /// In-place-add operator.
                virtual DenseElementIterator<DataType_> & operator+= (const unsigned long step)
                {
                    CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                    _index += step;

                    return *this;
                }

                /// Dereference operator that returns an assignable reference.
                virtual DataType_ & operator* ()
                {
                    CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                    return _vector._imp->elements[_vector._imp->stepsize * _index + _vector._imp->offset];
                }

                /// Dereference operator that returns an unassignable reference.
                virtual const DataType_ & operator* () const
                {
                    CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                    return _vector._imp->elements[_vector._imp->stepsize * _index + _vector._imp->offset];
                }

                /// Less-than operator.
                virtual bool operator< (const IteratorBase<DataType_, Vector<DataType_> > & other) const
                {
                    return _index < other.index();
                }

                /// Equality operator.
                virtual bool operator== (const IteratorBase<DataType_, Vector<DataType_> > & other) const
                {
                    return ((&_vector == other.parent()) && (_index == other.index()));
                }

                /// Inequality operator.
                virtual bool operator!= (const IteratorBase<DataType_, Vector<DataType_> > & other) const
                {
                    return ((&_vector != other.parent()) || (_index != other.index()));
                }

                /// \}

                /// \name IteratorTraits interface
                /// \{

                /// Returns our index.
                virtual unsigned long index() const
                {
                    return _index;
                }

                /// Returns a pointer to our parent container.
                virtual const Vector<DataType_> * parent() const
                {
                    return &_vector;
                }

                /// \}
        };
}

#endif
