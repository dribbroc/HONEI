/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_BANDED_MATRIX_IMPL_HH
#define LIBLA_GUARD_BANDED_MATRIX_IMPL_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_vector-impl.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/assertion.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <set>

namespace honei
{
    // Forward declarations
    template <typename DataType_> struct BandIteratorImplementation;
    template <typename DataType_> struct ConstBandIteratorImplementation;
    template <typename DataType_> struct ConstNonZeroBandIteratorImplementation;
    template <typename DataType_> struct NonZeroBandIteratorImplementation;

    template <typename DataType_> struct Implementation<BandedMatrix<DataType_> >
    {
        /// Array of pointers to our band-data.
        SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > bands;

        /// Our size.
        unsigned long size;

        /// Our zero element.
        static const DataType_ zero_element;

        /// Our zero vector.
        DenseVector<DataType_> zero_vector;

        /// Our set of non-zero bands
        std::set<unsigned long> non_zero_bands;

        Implementation(unsigned long size) :
            bands(2 * size - 1),
            size(size),
            zero_vector(size, DataType_(0))
        {
            non_zero_bands.insert(2 * size - 1);
        }
    };

    template <typename DataType_>
    BandedMatrix<DataType_>::BandedMatrix(unsigned long size) :
        PrivateImplementationPattern<BandedMatrix<DataType_>, Shared>(new Implementation<BandedMatrix<DataType_> >(size))
    {
        CONTEXT("When creating BandedMatrix:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    BandedMatrix<DataType_>::BandedMatrix(unsigned long size, const DenseVector<DataType_> & diagonal) :
        PrivateImplementationPattern<BandedMatrix<DataType_>, Shared>(new Implementation<BandedMatrix<DataType_> >(size))
    {
        CONTEXT("When creating BandedMatrix with initial band:");
        ASSERT(size > 0, "size is zero!");

        if (diagonal.size() != size)
            throw VectorSizeDoesNotMatch(diagonal.size(), size);

        this->_imp->bands[size - 1].reset(new DenseVector<DataType_>(diagonal));
        this->_imp->non_zero_bands.insert(size - 1);
    }

    template <typename DataType_>
    BandedMatrix<DataType_>::BandedMatrix(const BandedMatrixQ1<DataType_> & src) :
        PrivateImplementationPattern<BandedMatrix<DataType_>, Shared>(new Implementation<BandedMatrix<DataType_> >(src.size()))
    {
            this->insert_band(0, src.band(DD));
            this->insert_band(1, src.band(DU));
            this->insert_band(src.root() - 1, src.band(UL));
            this->insert_band(src.root(), src.band(UD));
            this->insert_band(src.root() + 1, src.band(UU));
            this->insert_band(-1, src.band(DL));
            this->insert_band(-src.root() - 1, src.band(LL));
            this->insert_band(-src.root(), src.band(LD));
            this->insert_band(-src.root() + 1, src.band(LU));
    }

    template <typename DataType_>
    BandedMatrix<DataType_>::~BandedMatrix()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_>
    BandedMatrix<DataType_>::begin_elements()
    {
        return ElementIterator(*this, 0);
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_>
    BandedMatrix<DataType_>::end_elements()
    {
        return ElementIterator(*this, this->_imp->size * this->_imp->size);
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>
    BandedMatrix<DataType_>::begin_elements() const
    {
        return ConstElementIterator(*this, 0);
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>
    BandedMatrix<DataType_>::end_elements() const
    {
        return ConstElementIterator(*this, this->_imp->size * this->_imp->size);
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::begin_bands()
    {
        return BandIterator(new BandIteratorImplementation<DataType_>(this->_imp, 0));
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::band_at(unsigned long index)
    {
        return BandIterator(new BandIteratorImplementation<DataType_>(this->_imp, index));
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::end_bands()
    {
        return BandIterator(new BandIteratorImplementation<DataType_>(this->_imp, 2 * this->_imp->size - 1));
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::begin_bands() const
    {
        return ConstBandIterator(new ConstBandIteratorImplementation<DataType_>(this->_imp, 0));
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::band_at(unsigned long index) const
    {
        return ConstBandIterator(new ConstBandIteratorImplementation<DataType_>(this->_imp, index));
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::end_bands() const
    {
        return ConstBandIterator(new ConstBandIteratorImplementation<DataType_>(this->_imp, 2 * this->_imp->size - 1));
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::begin_non_zero_bands()
    {
        return BandIterator(new NonZeroBandIteratorImplementation<DataType_>(this->_imp, this->_imp->non_zero_bands.begin()));
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::end_non_zero_bands()
    {
        return BandIterator(new NonZeroBandIteratorImplementation<DataType_>(this->_imp, this->_imp->non_zero_bands.find(2 * this->_imp->size - 1)));
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::begin_non_zero_bands() const
    {
        return ConstBandIterator(new ConstNonZeroBandIteratorImplementation<DataType_>(this->_imp, this->_imp->non_zero_bands.begin()));
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>
    BandedMatrix<DataType_>::end_non_zero_bands() const
    {
        return ConstBandIterator(new ConstNonZeroBandIteratorImplementation<DataType_>(this->_imp, this->_imp->non_zero_bands.find(2 * this->_imp->size - 1)));
    }

    template <typename DataType_>
    unsigned long
    BandedMatrix<DataType_>::columns() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long
    BandedMatrix<DataType_>::rows() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long
    BandedMatrix<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long
    BandedMatrix<DataType_>::non_zero_bands() const
    {
        return this->_imp->non_zero_bands.size();
    }

    template <typename DataType_>
    void
    BandedMatrix<DataType_>::insert_band(signed long index, const DenseVector<DataType_> & vector)
    {
        if (this->_imp->size != vector.size())
        {
            throw VectorSizeDoesNotMatch(this->_imp->size, vector.size());
        }

        if (! this->_imp->bands[index + this->_imp->size - 1])
        {
            this->_imp->non_zero_bands.insert(index + this->_imp->size - 1);
        }
        std::tr1::shared_ptr<DenseVector<DataType_> > temp(new DenseVector<DataType_>(vector));
        this->_imp->bands[index + this->_imp->size - 1] = temp;

    }

    template <typename DataType_>
    DenseVector<DataType_> &
    BandedMatrix<DataType_>::band(signed long index) const
    {
        CONTEXT("When retrieving band '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");
        ASSERT(std::abs(index) < (signed)this->_imp->size, "index out of bounds!");

        if (! this->_imp->bands[index + this->_imp->size - 1])
        {
            this->_imp->bands[index + this->_imp->size - 1].reset(new DenseVector<DataType_>(this->_imp->size, DataType_(0)));
            this->_imp->non_zero_bands.insert(index + this->_imp->size - 1);
        }

        return *this->_imp->bands[index + this->_imp->size - 1];
    }

    template <typename DataType_>
    DenseVector<DataType_> &
    BandedMatrix<DataType_>::band_unsigned(unsigned long index)
    {
        CONTEXT("When retrieving unsigned band '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");
        ASSERT(index < 2 * this->_imp->size - 1, "index '" + stringify(index) + "' out of bounds!");

        if (! this->_imp->bands[index])
        {
            this->_imp->bands[index].reset(new DenseVector<DataType_>(this->_imp->size, DataType_(0)));
            this->_imp->non_zero_bands.insert(index);
        }

        return *this->_imp->bands[index];
    }

    template <typename DataType_>
    void BandedMatrix<DataType_>::lock(LockMode mode) const
    {
        for (unsigned long i(0) ; i < 2 * this->_imp->size - 1 ; ++i)
        {
            if (this->_imp->bands[i])
            {
                this->_imp->bands[i]->lock(mode);
            }
        }
    }

    template <typename DataType_>
    void BandedMatrix<DataType_>::unlock(LockMode mode) const
    {
        for (unsigned long i(0) ; i < 2 * this->_imp->size - 1 ; ++i)
        {
            if (this->_imp->bands[i])
            {
                this->_imp->bands[i]->unlock(mode);
            }
        }
    }

    template <typename DataType_>
    BandedMatrix<DataType_>
    BandedMatrix<DataType_>::copy() const
    {
        CONTEXT("When creating copy() of a BandedMatrix:");
        BandedMatrix result(this->_imp->size);

        for (unsigned long i(0) ; i < 2 * this->_imp->size - 1 ; ++i)
        {
            if (this->_imp->bands[i])
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > temp(new DenseVector<DataType_>(
                            this->_imp->bands[i]->copy()));
                result._imp->bands[i] = temp;
                result._imp->non_zero_bands.insert(i);
            }
        }
        std::copy(this->_imp->non_zero_bands.begin(), this->_imp->non_zero_bands.end(),
                std::inserter(result._imp->non_zero_bands, result._imp->non_zero_bands.begin()));

        return result;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Banded, container::Matrix, DataType_> >
    {
        std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > matrix;

        unsigned long index;

        Implementation(const BandedMatrix<DataType_> & matrix, unsigned long index) :
            matrix(matrix._imp),
            index(index)
        {
        }

        /// Returns true if we're below the diagonal.
        inline bool is_lower() const
        {
            return index / matrix->size > index % matrix->size;
        }

        /// Returns the band-index for the current element.
        inline signed long band_index() const
        {
            signed long result(this->index % (this->matrix->size + 1));

            return is_lower() ? result - (this->matrix->size + 1) : result;
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_>::ElementIterator(const BandedMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::Banded, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Banded, container::Matrix, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_>::ElementIterator(const ElementIterator<storage::Banded, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ElementIterator<storage::Banded, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Banded, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_> &
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator= (
            const ElementIterator<storage::Banded, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->matrix = other._imp->matrix;
        this->_imp->index = other._imp->index;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_> &
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing iterator by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Banded, container::Matrix, DataType_> &
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->index) + "':");

        if (! this->_imp->matrix->bands[this->_imp->band_index() + this->_imp->matrix->size - 1])
        {
            this->_imp->matrix->bands[this->_imp->band_index() + this->_imp->matrix->size - 1].reset(new DenseVector<DataType_>(this->_imp->matrix->size, DataType_(0)));
            this->_imp->matrix->non_zero_bands.insert(this->_imp->band_index() + this->_imp->matrix->size - 1);
        }

        return (*this->_imp->matrix->bands[this->_imp->band_index() + this->_imp->matrix->size - 1])[row()];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator< (const ElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator== (const ElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() == other._imp->matrix.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Banded, container::Matrix, DataType_>::operator!= (const ElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() != other._imp->matrix.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Banded, container::Matrix, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Banded, container::Matrix, DataType_>::column() const
    {
        return this->_imp->index % this->_imp->matrix->size;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Banded, container::Matrix, DataType_>::row() const
    {
        return this->_imp->index / this->_imp->matrix->size;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Banded, container::Matrix, DataType_> >
    {
        std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > matrix;

        unsigned long index;

        Implementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix, unsigned long index) :
            matrix(matrix),
            index(index)
        {
        }

        /// Returns true if we're below the diagonal.
        inline bool is_lower() const
        {
            return index / matrix->size > index % matrix->size;
        }

        /// Returns the band-index for the current element.
        inline signed long band_index() const
        {
            signed long result(index % (matrix->size + 1));

            return is_lower() ? result - (matrix->size + 1) : result;
        }

    };

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::ConstElementIterator(const BandedMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Banded, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Banded, container::Matrix, DataType_> >(matrix._imp, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::ConstElementIterator(const ConstElementIterator<storage::Banded, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Banded, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Banded, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::ConstElementIterator(const ElementIterator<storage::Banded, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Banded, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Banded, container::Matrix, DataType_> >(other._imp->matrix, other._imp->index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_> &
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator= (
            const ConstElementIterator<storage::Banded, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->matrix = other._imp->matrix;
        this->_imp->index = other._imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_> &
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing iterator by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Banded, container::Matrix, DataType_> &
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->index) + "':");

        if (! this->_imp->matrix->bands[this->_imp->band_index() + this->_imp->matrix->size - 1])
        {
            return Implementation<BandedMatrix<DataType_> >::zero_element;
        }
        else
        {
            return (*this->_imp->matrix->bands[this->_imp->band_index() + this->_imp->matrix->size - 1])[row()];
        }
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator< (
            const ConstElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator== (
            const ConstElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() == other._imp->matrix.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::operator!= (
            const ConstElementIterator<storage::Banded, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() != other._imp->matrix.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::column() const
    {
        return this->_imp->index % this->_imp->matrix->size;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Banded, container::Matrix, DataType_>::row() const
    {
        return this->_imp->index / this->_imp->matrix->size;
    }

    template <typename DataType_> struct Implementation<BandIterator<type::Banded, DataType_> >
    {
        std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > matrix;

        unsigned long index;

        DenseVector<DataType_> * band;

        Implementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix, unsigned long index) :
            matrix(matrix),
            index(index),
            band(0)
        {
        }

        virtual ~Implementation()
        {
        }

        virtual void advance() = 0;

        virtual Implementation * clone() const = 0;

        virtual bool exists() const = 0;

        virtual Implementation<ConstBandIterator<type::Banded, DataType_> > * make_const() const = 0;
    };

    template <typename DataType_> struct BandIteratorImplementation :
        public Implementation<BandIterator<type::Banded, DataType_> >
    {
        BandIteratorImplementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix, unsigned long index) :
            Implementation<BandIterator<type::Banded, DataType_> >(matrix, index)
        {
        }

        virtual ~BandIteratorImplementation()
        {
        }

        virtual void advance()
        {
            ++this->index;
        }

        virtual BandIteratorImplementation * clone() const
        {
            return new BandIteratorImplementation(*this);
        }

        virtual bool exists() const
        {
            return this->matrix->bands[this->index];
        }

        virtual Implementation<ConstBandIterator<type::Banded, DataType_> > * make_const() const
        {
            return new ConstBandIteratorImplementation<DataType_>(this->matrix, this->index);
        }
    };

    template <typename DataType_> struct NonZeroBandIteratorImplementation :
        public Implementation<BandIterator<type::Banded, DataType_> >
    {
        std::set<unsigned long>::const_iterator iterator;

        NonZeroBandIteratorImplementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix,
                const std::set<unsigned long>::const_iterator & iterator) :
            Implementation<BandIterator<type::Banded, DataType_> >(matrix, *iterator),
            iterator(iterator)
        {
            CONTEXT("When creating NonZeroBandIteratorImplementation:");
        }

        virtual ~NonZeroBandIteratorImplementation()
        {
        }

        virtual void advance()
        {
            CONTEXT("When advancing NonZeroBandIteratorImplementation:");

            ++this->iterator;
            this->index = *this->iterator;
        }

        virtual NonZeroBandIteratorImplementation * clone() const
        {
            return new NonZeroBandIteratorImplementation(*this);
        }

        virtual bool exists() const
        {
            return true;
        }

        virtual Implementation<ConstBandIterator<type::Banded, DataType_> > * make_const() const
        {
            return new ConstNonZeroBandIteratorImplementation<DataType_>(this->matrix, this->iterator);
        }
    };

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>::BandIterator(Implementation<BandIterator<type::Banded, DataType_> > * imp) :
        PrivateImplementationPattern<BandIterator<type::Banded, DataType_>, Single>(imp)
    {
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>::BandIterator(const BandIterator<type::Banded, DataType_> & other) :
        PrivateImplementationPattern<BandIterator<type::Banded, DataType_>, Single>(other._imp->clone())
    {
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_>::~BandIterator()
    {
    }

    template <typename DataType_>
    BandIterator<type::Banded, DataType_> &
    BandIterator<type::Banded, DataType_>::operator++ ()
    {
        this->_imp->advance();

        return *this;
    }

    template <typename DataType_>
    bool
    BandIterator<type::Banded, DataType_>::operator== (const BandIterator<type::Banded, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() == other._imp->matrix.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    BandIterator<type::Banded, DataType_>::operator!= (const BandIterator<type::Banded, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() != other._imp->matrix.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    DenseVector<DataType_>
    BandIterator<type::Banded, DataType_>::operator* () const
    {
        if (!this->_imp->matrix->bands[this->_imp->index])
        {
            this->_imp->matrix->bands[this->_imp->index].reset(new DenseVector<DataType_>(this->_imp->matrix->size, DataType_(0)));
            this->_imp->matrix->non_zero_bands.insert(this->_imp->index);
        }

        return *this->_imp->matrix->bands[this->_imp->index];
    }

    template <typename DataType_>
    bool
    BandIterator<type::Banded, DataType_>::exists() const
    {
        return this->_imp->exists();
    }

    template <typename DataType_>
    unsigned long
    BandIterator<type::Banded, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_> struct Implementation<ConstBandIterator<type::Banded, DataType_> >
    {
        std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > matrix;

        unsigned long index;

        Implementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix, unsigned long index) :
            matrix(matrix),
            index(index)
        {
        }

        virtual ~Implementation()
        {
        }

        virtual void advance() = 0;

        virtual Implementation * clone() const = 0;

        virtual bool exists() = 0;

    };

    template <typename DataType_> struct ConstBandIteratorImplementation :
        public Implementation<ConstBandIterator<type::Banded, DataType_> >
    {
        ConstBandIteratorImplementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix, unsigned long index) :
            Implementation<ConstBandIterator<type::Banded, DataType_> >(matrix, index)
        {
        }

        virtual ~ConstBandIteratorImplementation()
        {
        }

        virtual void advance()
        {
            ++this->index;
        }

        virtual ConstBandIteratorImplementation * clone() const
        {
            return new ConstBandIteratorImplementation(*this);
        }

        virtual bool exists()
        {
            return this->matrix->bands[this->index];
        }
    };

    template <typename DataType_> struct ConstNonZeroBandIteratorImplementation :
        public Implementation<ConstBandIterator<type::Banded, DataType_> >
    {
        std::set<unsigned long>::const_iterator iterator;

        ConstNonZeroBandIteratorImplementation(const std::tr1::shared_ptr<Implementation<BandedMatrix<DataType_> > > & matrix,
                const std::set<unsigned long>::const_iterator & iterator) :
            Implementation<ConstBandIterator<type::Banded, DataType_> >(matrix, *iterator),
            iterator(iterator)
        {
            CONTEXT("When creating ConstNonZeroBandIteratorImplementation:");
        }

        virtual ~ConstNonZeroBandIteratorImplementation()
        {
        }

        virtual void advance()
        {
            CONTEXT("When advancing NonZeroBandIteratorImplementation:");
            ASSERT(this->iterator != this->matrix->non_zero_bands.end(), "Incrementing end iterator!");

            ++this->iterator;
            this->index = *this->iterator;
        }

        virtual ConstNonZeroBandIteratorImplementation * clone() const
        {
            return new ConstNonZeroBandIteratorImplementation(*this);
        }

        virtual bool exists()
        {
            return true;
        }
    };

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>::ConstBandIterator(Implementation<ConstBandIterator<type::Banded, DataType_> > * imp) :
        PrivateImplementationPattern<ConstBandIterator<type::Banded, DataType_>, Single>(imp)
    {
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>::ConstBandIterator(const ConstBandIterator<type::Banded, DataType_> & other) :
        PrivateImplementationPattern<ConstBandIterator<type::Banded, DataType_>, Single>(other._imp->clone())
    {
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>::ConstBandIterator(const BandIterator<type::Banded, DataType_> & other) :
        PrivateImplementationPattern<ConstBandIterator<type::Banded, DataType_>, Single>(other._imp->make_const())
    {
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_>::~ConstBandIterator()
    {
    }

    template <typename DataType_>
    ConstBandIterator<type::Banded, DataType_> &
    ConstBandIterator<type::Banded, DataType_>::operator++ ()
    {
        this->_imp->advance();

        return *this;
    }

    template <typename DataType_>
    bool
    ConstBandIterator<type::Banded, DataType_>::operator== (const ConstBandIterator<type::Banded, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() == other._imp->matrix.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ConstBandIterator<type::Banded, DataType_>::operator!= (const ConstBandIterator<type::Banded, DataType_> & other) const
    {
        return ((this->_imp->matrix.get() != other._imp->matrix.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    DenseVector<DataType_>
    ConstBandIterator<type::Banded, DataType_>::operator* () const
    {
        if (! this->_imp->matrix->bands[this->_imp->index])
        {
            return this->_imp->matrix->zero_vector;
        }
        else
        {
            return *this->_imp->matrix->bands[this->_imp->index];
        }
    }

    template <typename DataType_>
    DenseVector<DataType_>
    ConstBandIterator<type::Banded, DataType_>::operator-> () const
    {
        return *(*this);
    }

    template <typename DataType_>
    bool
    ConstBandIterator<type::Banded, DataType_>::exists() const
    {
        return this->_imp->exists();
    }

    template <typename DataType_>
    unsigned long
    ConstBandIterator<type::Banded, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    bool
    operator== (const BandedMatrix<DataType_> & a, const BandedMatrix<DataType_> & b)
    {
        if (a.size() != b.size())
            throw MatrixSizeDoesNotMatch(a.size(), b.size());

        for (typename BandedMatrix<DataType_>::ConstBandIterator i(a.begin_bands()), i_end(a.end_bands()),
                j(b.begin_bands()) ; i != i_end ; ++i, ++j)
        {
            for (typename DenseVector<DataType_>::ConstElementIterator x(i->begin_elements()), x_end(i->end_elements()),
                    y(j->begin_elements()) ; x < x_end ; ++x, ++y)
            {
                if (std::fabs(*x - *y) > std::numeric_limits<DataType_>::epsilon())
                    return false;
            }
        }

        return true;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const BandedMatrix<DataType_> & b)
    {
        lhs << "[";
        for (typename BandedMatrix<DataType_>::ConstBandIterator i(b.begin_bands()), i_end(b.end_bands()) ;
                i != i_end ; ++i)
        {
            lhs << "  " << *i;
        }
        lhs << "]" << std::endl;

        return lhs;
    }
}

#endif
