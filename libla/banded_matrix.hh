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
#include <libla/vector_iterator.hh>
#include <libutil/shared_array.hh>
#include <libutil/log.hh>

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
            /// Array of pointers to our bands.
            SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > _bands;

            /// Our size.
            unsigned long _size;

            /// Our zero element.
            static const DataType_ _zero_element;

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class BandedElementIterator;

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class BandedElementIterator<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Matrix<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename MutableMatrix<DataType_>::ElementIterator ElementIterator;

            /// Type of the const iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, const DataType_> ConstVectorIterator;

            /// Type of the iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, DataType_> VectorIterator;

            /// Constructors
            /// \{

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

            BandedMatrix(unsigned long size, DenseVector<DataType_> * diagonal) :
                _bands(new std::tr1::shared_ptr<DenseVector<DataType_> >[2 * size + 1]),
                _size(size)
            {
                if (diagonal->size() != size)
                    throw VectorSizeDoesNotMatch(diagonal->size(), size);

                _bands[0].reset(diagonal);
            }

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new BandedElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new BandedElementIterator<DataType_>(*this, _size * _size));
            }

            /// Returns the number of our columns.
            virtual unsigned long columns() const
            {
                return _size;
            }

            /// Returns the number of our rows.
            virtual unsigned long rows() const
            {
                return _size;
            }

            /// Returns a band-vector by index.
            DenseVector<DataType_> & band(unsigned long index) const
            {
                if (! _bands[index])
                    _bands[index].reset(new DenseVector<DataType_>(index < _size ? _size - index : index + 1 - _size,
                                static_cast<DataType_>(0)));
                return *_bands[index];
            }

#if 0
            /// Returns iterator pointing to the first non-zero band of the matrix.
            ConstVectorIterator begin_bands() const
            {
                return ConstVectorIterator(new VectorIteratorImpl<DataType_>(*this, 0));
            }

            /// Return iterator pointing behind the last non-zero band of the matrix.
            ConstVectorIterator end_bands() const
            {
                return ConstVectorIterator(new VectorIteratorImpl<DataType_>(*this, 2 * _size + 1));
            }
#endif
    };

    template <typename DataType_> const DataType_ BandedMatrix<DataType_>::_zero_element(0);

    /**
     * BandedMatrix::BandedElementIterator is a simple iterator implementation for banded matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class BandedMatrix<DataType_>::BandedElementIterator<DataType_> :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
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
            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            BandedElementIterator(const BandedMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            BandedElementIterator(BandedElementIterator<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// \}

            /// Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual MatrixElementIterator & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Dereference operator thar returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                if (! _matrix._bands[_band_index()])
                {
                    _matrix._bands[_band_index()].reset(new DenseVector<DataType_>(_matrix._size, static_cast<DataType_>(0)));
                }

                return (*_matrix._bands[_band_index()])[(_lower() ? _index : _index - _band_class()) % _matrix._size];
            }

            /// Dereference operator thar returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                if (! _matrix._bands[_band_index()])
                    return _matrix._zero_element;
                else
                    return (*_matrix._bands[_band_index()])[(_lower() ? _index : _index - _band_class()) % _matrix._size];
            }

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_matrix == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// \}

            /// IteratorTraits interface
            /// \{

            /// Returns our index.
            virtual unsigned long index() const
            {
                return _index;
            }

            /// Returns our column index.
            virtual unsigned long column() const
            {
                return _index % _matrix._size;
            }

            /// Returns our row index.
            virtual unsigned long row() const
            {
                return _index / _matrix._size;
            }

            /// Returns a pointer to our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }

            /// \}
    };

#if 0
    /**
     * A BandedMatrix::VectorIteratorImpl is a simple vector iterator implementation for banded matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class BandedMatrix<DataType_>::VectorIteratorImpl<const DataType_> :
        public VectorIteratorImplBase<DataType_, const Vector<DataType_> >
    {
        private:
            /// Our matrix.
            const BandedMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index of the matrix's band.
             **/
            VectorIteratorImpl(const BandedMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            VectorIteratorImpl(VectorIteratorImpl<const DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual VectorIteratorImpl<const DataType_> & operator++ ()
            {
                ++_index;
                for ( ; (_index <= 2 * _matrix._size + 1) && (! _matrix._bands[_index]) ; ++_index)
                    ;

                return *this;
            }

            /// Postincrement operator.
            virtual VectorIteratorImpl<const DataType_> operator++ (int)
            {
                VectorIteratorImpl<const DataType_> result(*this);
                ++_index;
                for ( ; (_index <= 2 * _matrix._size + 1) && (! _matrix._bands[_index]) ; ++_index)
                    ;

                return result;
            }

            /// Equality operator.
            virtual bool operator== (const VectorIteratorImplBase<DataType_, const Vector<DataType_> > & other) const
            {
                return (&_matrix == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const VectorIteratorImplBase<DataType_, const Vector<DataType_> > & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual const Vector<DataType_> & operator* () const
            {
                return *_matrix._bands[_index];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };
#endif
}

#endif
