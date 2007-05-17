/* vim: set sw=4 sts=4 et nofoldenable : */

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

#ifndef LIBLA_GUARD_BANDED_MATRIX_HH
#define LIBLA_GUARD_BANDED_MATRIX_HH 1

#include <libla/element_iterator.hh>
#include <libla/matrix.hh>
#include <libla/dense_vector.hh>
#include <libutil/shared_array.hh>

#include <string.h>
#include <iterator>
#include <tr1/memory>

/**
 * \file
 *
 * Implementation of BandedMatrix and related classes.
 *
 * \ingroup grpmatrix
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A BandedMatrix is a square matrix with O(size) non-zero bands which keeps its data
     * non-sequential.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class BandedMatrix :
        public Matrix<DataType_>
    {
        private:
            /// Pointer to our bands.
            SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > _bands;

            /// Our size.
            unsigned long _size;

            /// Our zero element.
            static DataType_ _zero_element;

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_> ElementIterator;

            /// Type of the iterator over our row/column/band/diagonal vectors.
//            typedef VectorIteratorBase<Matrix<DataType_>, Vector<DataType_> > VectorIterator;

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             **/
            BandedMatrix(unsigned long size) :
                _bands(new std::tr1::shared_ptr<DenseVector<DataType_> >[2 * size + 1]),
                _size(size)
            {
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, _size * _size));
            }

            /// Returns our columns.
            virtual unsigned long columns() const
            {
                return _size;
            }

            /// Returns our rows.
            virtual unsigned long rows() const
            {
                return _size;
            }

            /// Return a band by index.
            DenseVector<DataType_> & band(unsigned long index) const
            {
                if (! _bands[index])
                    _bands[index].reset(new DenseVector<DataType_>(index < _size ? _size - index : index + 1 - _size));
                return *_bands[index];
            }
    };

    template <typename DataType_> DataType_ BandedMatrix<DataType_>::_zero_element(0);

    /**
     * A BandedMatrix::ElementIteratorImpl is a simple iterator implementation for banded matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class BandedMatrix<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Matrix<DataType_>, DataType_>
    {
        private:
            /// Our matrix.
            const BandedMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

            /// Returns the band-class for the current element.
            inline unsigned long _band_class() const
            {
                return _index % (_matrix._size + 1);
            }

            /// Returns true if we're below the diagonal.
            inline bool _lower() const
            {
                return _index % _matrix._size < _band_class();
            }

            /// Returns the band-index for the current element.
            inline unsigned long _band_index() const
            {
                return (_lower() ? _matrix._size - 2 + _band_class() : _band_class());
            }

            /// Return the band-size for the current element.
            inline unsigned long _band_size() const
            {
                return (_band_index() < _matrix._size ?
                        _matrix._size - _band_index() : _band_index() + 1 - _matrix._size);
            }

        public:
            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            ElementIteratorImpl(const BandedMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            ElementIteratorImpl(ElementIteratorImpl<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<DataType_> operator++ (int)
            {
                ElementIteratorImpl<DataType_> result(*this);
                ++_index;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Matrix<DataType_>, DataType_> & other) const
            {
                return (&_matrix == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Matrix<DataType_>, DataType_> & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                if (! _matrix._bands[_band_index()])
                    _matrix._bands[_band_index()].reset(new DenseVector<DataType_>(_band_size(), static_cast<DataType_>(0)));

                return (*_matrix._bands[_band_index()])[(_lower() ? _index : _index - _band_class()) % _matrix._size];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns our column.
            virtual const unsigned long column() const
            {
                return _index % _matrix._size;
            }

            /// Returns our row.
            virtual const unsigned long row() const
            {
                return _index / _matrix._size;
            }

            /// Returns our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };
}

#endif
