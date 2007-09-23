/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 *
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
#include <libla/vector_error.hh>
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
namespace honei
{
    /**
     * \brief BandedMatrix is a square matrix with O(size) non-zero bands which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class BandedMatrix :
        public Matrix<DataType_>
    {
        private:
            /// Array of pointers to our band-data.
            SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > _bands;

            /// Our size.
            unsigned long _size;

            /// Our zero element.
            static const DataType_ _zero_element;

            /// Our zero vector.
            DenseVector<DataType_> _zero_vector;

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class BandedElementIterator;

            /// Our implementation of VectorIteratorBase.
            template <typename ElementType_> class BandIterator;

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class BandedElementIterator<DataType_>;
            friend class BandIterator<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Matrix<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename MutableMatrix<DataType_>::ElementIterator ElementIterator;

            /// Type of the iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, DenseVector<DataType_> > VectorIterator;

            /// Type of the const iterator over our vectors.
            typedef ConstVectorIteratorWrapper<DataType_, DenseVector<DataType_> > ConstVectorIterator;

            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             **/
            BandedMatrix(unsigned long size) :
                _bands(2 * size - 1),
                _size(size),
                _zero_vector(size, DataType_(0))
            {
                CONTEXT("When creating BandedMatrix:");
            }

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             * \param diagonal Diagonal of the new banded matrix.
             **/
            BandedMatrix(unsigned long size, const DenseVector<DataType_> & diagonal) :
                _bands(2 * size - 1),
                _size(size),
                _zero_vector(size, DataType_(0))
            {
                CONTEXT("When creating BandedMatrix with initial band:");
                if (diagonal.size() != size)
                    throw VectorSizeDoesNotMatch(diagonal.size(), size);

                _bands[size - 1].reset(new DenseVector<DataType_>(diagonal));
            }
            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new BandedElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new BandedElementIterator<DataType_>(*this, _size * _size));
            }

            /// Returns const iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new BandedElementIterator<DataType_>(*this, 0));
            }

            /// Returns const iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new BandedElementIterator<DataType_>(*this, _size * _size));
            }

            /// Returns iterator pointing to the first band of the matrix.
            VectorIterator begin_bands()
            {
                return VectorIterator(new BandIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing to a given band of the matrix.
            VectorIterator band_at(unsigned long index)
            {
                return VectorIterator(new BandIterator<DataType_>(*this, index));
            }

            /// Returns iterator pointing behind the last band of the matrix.
            VectorIterator end_bands()
            {
                return VectorIterator(new BandIterator<DataType_>(*this, 2 * _size - 1));
            }

            /// Returns iterator pointing to the first band of the matrix.
            ConstVectorIterator begin_bands() const
            {
                return ConstVectorIterator(new BandIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing to a given band of the matrix.
            ConstVectorIterator band_at(unsigned long index) const
            {
                return VectorIterator(new BandIterator<DataType_>(*this, index));
            }

            /// Returns iterator pointing behind the last band of the matrix.
            ConstVectorIterator end_bands() const
            {
                return ConstVectorIterator(new BandIterator<DataType_>(*this, 2 * _size - 1));
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

            /// Returns our size, equal to rows and columns.
            virtual unsigned long size() const
            {
                return _size;
            }

            /// Inserts a new Band in the matrix.
            void insert_band(signed long index, const DenseVector<DataType_> & vector)
            {
                if (_size != vector.size())
                {
                    throw VectorSizeDoesNotMatch(_size, vector.size());
                }

                std::tr1::shared_ptr<DenseVector<DataType_> > temp(new DenseVector<DataType_>(vector));
                _bands[index + _size - 1] = temp;

            }

            /// Returns a band-vector by index.
            DenseVector<DataType_> & band(signed long index) const
            {
                CONTEXT("When retrieving band '" + stringify(index) + "' of matrix of size '"
                        + stringify(_size) + "':");
                ASSERT(std::abs(index) < _size, "index out of bounds!");

                if (! _bands[index + _size - 1])
                    _bands[index + _size - 1].reset(new DenseVector<DataType_>(_size, DataType_(0)));

                return *_bands[index + _size - 1];
            }

            /// Returns a copy of the matrix.
            virtual BandedMatrix * copy() const
            {
                CONTEXT("When creating copy() of a BandedMatrix:");
                BandedMatrix * result(new BandedMatrix(_size));

                for (unsigned long i(0) ; i < 2 * _size - 1 ; ++i)
                {
                    if (_bands[i])
                    {
                        std::tr1::shared_ptr<DenseVector<DataType_> > temp(new DenseVector<DataType_>(
                                    _bands[i]->copy()));
                        result->_bands[i] = temp;
                    }
                }

                return result;
            }
    };

    template <typename DataType_> const DataType_ BandedMatrix<DataType_>::_zero_element(0);

    /**
     * \brief BandedMatrix::BandedElementIterator is a simple iterator implementation for banded matrices.
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

            /// Returns true if we're below the diagonal.
            inline bool _lower() const
            {
                return row() > column();
            }

            /// Returns the band-index for the current element.
            inline signed long _band_index() const
            {
                signed long result(_index % (_matrix._size + 1));

                return _lower() ? result - (_matrix._size + 1) : result;
            }

        public:
            /// \name Constructors
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

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual MatrixElementIterator & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

                ++_index;

                return *this;
            }

            /// In-place-add operator.
            virtual MatrixElementIterator & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                if (! _matrix._bands[_band_index() + _matrix._size - 1])
                {
                    _matrix._bands[_band_index() + _matrix._size - 1].reset(new DenseVector<DataType_>(_matrix._size, DataType_(0)));
                }

                return (*_matrix._bands[_band_index() + _matrix._size - 1])[row()];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                if (! _matrix._bands[_band_index() + _matrix._size - 1])
                {
                    return _matrix._zero_element;
                }
                else
                {
                    return (*_matrix._bands[_band_index() + _matrix._size - 1])[row()];
                }
            }

            /// Comparison operator for less-than.
            virtual bool operator< (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return _index < other.index();
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

            /// \name IteratorTraits interface
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

    template <> template <typename DataType_> class BandedMatrix<DataType_>::BandIterator<DataType_> :
            public VectorIteratorBase<DataType_, DenseVector<DataType_> >
    {
            private:
            /// Our parent matrix.
            const BandedMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/

            BandIterator(const BandedMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            BandIterator(BandIterator<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual VectorIteratorBase<DataType_, DenseVector<DataType_> > & operator++ ()
            {
                ++_index;

                return *this;
            }

            /// Equality operator.
            virtual bool operator== (const VectorIteratorBase<DataType_, DenseVector<DataType_> > & other) const
            {
                return ((&_matrix == other.parent()) && (_index == other.index()));
            }

            /// Inequality operator.
            virtual bool operator!= (const VectorIteratorBase<DataType_, DenseVector<DataType_> > & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator that returns an assignable reference.
            virtual DenseVector<DataType_> & operator* ()
            {
                if (!_matrix._bands[_index])
                {
                    _matrix._bands[_index].reset(new DenseVector<DataType_>(_matrix._size, DataType_(0)));
                }

                return (*_matrix._bands[_index]);
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DenseVector<DataType_> & operator* () const
            {
                if (! _matrix._bands[_index])
                {
                    return _matrix._zero_vector;
                }
                else
                {
                    return (*_matrix._bands[_index]);
                }
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns true if the referenced vector already exists.
            virtual bool exists() const
            {
                return _matrix._bands[_index];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns a pointer to our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }

            /// \}
    };


}

#endif
