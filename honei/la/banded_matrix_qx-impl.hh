/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBLA_GUARD_BANDED_MATRIX_QX_IMPL_HH
#define LIBLA_GUARD_BANDED_MATRIX_QX_IMPL_HH 1

#include <honei/la/banded_matrix_qx.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/assertion.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/la/band_type.hh>

#include <cmath>

namespace honei
{
    template <BandType BandType_, typename DataType_> struct Implementation<BandedMatrixQx<BandType_, DataType_> >
    {
    };

    template <typename DataType_> struct Implementation<BandedMatrixQx<Q1Type, DataType_> >
    {
        /// Array of pointers to our band-data.
        SharedArray<shared_ptr<DenseVector<DataType_> > > bands;

        /// Number of bands
        unsigned long band_count;

        /// Our size.
        unsigned long size;

        /// The square roote of our size.
        signed long root;

        /// Our zero element.
        static const DataType_ zero_element;

        /// Our zero vector.
        DenseVector<DataType_> zero_vector;

        Implementation(unsigned long size) :
            bands(9),
            band_count(9),
            size(size),
            root(static_cast<signed long>(std::sqrt(size))),
            zero_vector(size, DataType_(0))
        {
            unsigned long target;
            for (unsigned long level(0) ; level <= 50 ; ++level)
            {
                target = (unsigned long)pow((pow(2, level) + 1), 2);
                if (target >= size)
                    break;
            }
            ASSERT(size == target, "size of " + stringify(size) + " matches not ((2^L) +1)^2, L integer.");
        }

        Implementation(BandedMatrix<DataType_> & src) :
            bands(9),
            size(src.size()),
            root(static_cast<signed long>(std::sqrt(size))),
            zero_vector(size, DataType_(0))
        {
            ASSERT(src.size() > 0, "size is zero!");

            unsigned long target;
            for (unsigned long level(0) ; level <= 50 ; ++level)
            {
                target = (unsigned long)pow((pow(2, level) + 1), 2);
                if (target >= size)
                    break;
            }
            ASSERT(size == target, "size of " + stringify(size) + " matches not ((2^L) +1)^2, L integer.");


            bands[LL].reset(new DenseVector<DataType_>(src.band(-root - 1)));
            bands[LD].reset(new DenseVector<DataType_>(src.band(-root)));
            bands[LU].reset(new DenseVector<DataType_>(src.band(-root + 1)));
            bands[DL].reset(new DenseVector<DataType_>(src.band(-1)));
            bands[DD].reset(new DenseVector<DataType_>(src.band(0)));
            bands[DU].reset(new DenseVector<DataType_>(src.band(1)));
            bands[UL].reset(new DenseVector<DataType_>(src.band(root - 1)));
            bands[UD].reset(new DenseVector<DataType_>(src.band(root)));
            bands[UU].reset(new DenseVector<DataType_>(src.band(root + 1)));
        }

        Implementation(SparseMatrixELL<DataType_> & src) :
            bands(9),
            size(src.rows()),
            root(static_cast<signed long>(std::sqrt(size))),
            zero_vector(size, DataType_(0))
        {
            ASSERT(src.rows() > 0, "size is zero!");
            ASSERT(src.rows() == src.columns(), "Matrix is not squared!");

            unsigned long target;
            for (unsigned long level(0) ; level <= 50 ; ++level)
            {
                target = (unsigned long)pow((pow(2, level) + 1), 2);
                if (target >= size)
                    break;
            }
            ASSERT(size == target, "size of " + stringify(size) + " matches not ((2^L) +1)^2, L integer.");

            bands[LL].reset(new DenseVector<DataType_>(src.rows()));
            bands[LD].reset(new DenseVector<DataType_>(src.rows()));
            bands[LU].reset(new DenseVector<DataType_>(src.rows()));
            bands[DL].reset(new DenseVector<DataType_>(src.rows()));
            bands[DD].reset(new DenseVector<DataType_>(src.rows()));
            bands[DU].reset(new DenseVector<DataType_>(src.rows()));
            bands[UL].reset(new DenseVector<DataType_>(src.rows()));
            bands[UD].reset(new DenseVector<DataType_>(src.rows()));
            bands[UU].reset(new DenseVector<DataType_>(src.rows()));

            for (unsigned long i(0) ; i < src.rows() ; ++i)
            {
                if (i + root + 1 < src.rows())
                    (*bands[LL])[i + root + 1] = src(i + root + 1, i);

                if (i + root < src.rows())
                    (*bands[LD])[i + root] = src(i + root, i);

                if (i + root - 1 < src.rows())
                    (*bands[LU])[i + root - 1] = src(i + root - 1, i);

                if (i + 1 < src.rows())
                    (*bands[DL])[i + 1] = src(i + 1, i);

                (*bands[DD])[i] = src(i, i);

                if (i + 1 < src.rows())
                    (*bands[DU])[i] = src(i, i + 1);

                if (i + root - 1 < src.rows())
                    (*bands[UL])[i] = src(i, i + root - 1);

                if (i + root < src.rows())
                    (*bands[UD])[i] = src(i, i + root);

                if (i + root + 1 < src.rows())
                    (*bands[UU])[i] = src(i, i + root + 1);
            }
        }

        const DataType_ & operator() (unsigned long row, unsigned long column) const
        {
            long index(column - row);

            if (index == -root - 1)
                return (*bands[LL])[row];
            else if (index == -root)
                return (*bands[LD])[row];
            else if (index == -root + 1)
                return (*bands[LU])[row];
            else if (index == -1)
                return (*bands[DL])[row];
            else if (index == 0)
                return (*bands[DD])[row];
            else if (index == 1)
                return (*bands[DU])[row];
            else if (index == root - 1)
                return (*bands[UL])[row];
            else if (index == root)
                return (*bands[UD])[row];
            else if (index == root + 1)
                return (*bands[UU])[row];
            else
                return zero_element;
        }

        DenseVectorRange<DataType_> band_range(unsigned long index) const
        {
            switch(index)
            {
                case LL:
                    return DenseVectorRange<DataType_>(*bands[LL], size - root - 1, root + 1);
                    break;
                case LD:
                    return DenseVectorRange<DataType_>(*bands[LD], size - root, root);
                    break;
                case LU:
                    return DenseVectorRange<DataType_>(*bands[LU], size - root + 1, root - 1 );
                    break;
                case DL:
                    return DenseVectorRange<DataType_>(*bands[DL], size - 1, 1);
                    break;
                case DD:
                    return DenseVectorRange<DataType_>(*bands[DD], size, 0);
                    break;
                case DU:
                    return DenseVectorRange<DataType_>(*bands[DU], size - 1, 0);
                    break;
                case UL:
                    return DenseVectorRange<DataType_>(*bands[UL], size - root + 1, 0);
                    break;
                case UD:
                    return DenseVectorRange<DataType_>(*bands[UD], size - root, 0);
                    break;
                case UU:
                    return DenseVectorRange<DataType_>(*bands[UU], size - root - 1, 0);
                    break;
                default:
                    return DenseVectorRange<DataType_>(zero_vector, size, 0);
            }
        }

        void lock(LockMode mode) const
        {
            bands[LL]->lock(mode);
            bands[LD]->lock(mode);
            bands[LU]->lock(mode);
            bands[DL]->lock(mode);
            bands[DD]->lock(mode);
            bands[DU]->lock(mode);
            bands[UL]->lock(mode);
            bands[UD]->lock(mode);
            bands[UU]->lock(mode);
        }

        void unlock(LockMode mode) const
        {
            bands[LL]->unlock(mode);
            bands[LD]->unlock(mode);
            bands[LU]->unlock(mode);
            bands[DL]->unlock(mode);
            bands[DD]->unlock(mode);
            bands[DU]->unlock(mode);
            bands[UL]->unlock(mode);
            bands[UD]->unlock(mode);
            bands[UU]->unlock(mode);
        }

        BandedMatrixQx<Q1Type, DataType_> copy() const
        {
            BandedMatrixQx<Q1Type, DataType_> result(size,
                    (bands[0])->copy(),
                    (bands[1])->copy(),
                    (bands[2])->copy(),
                    (bands[3])->copy(),
                    (bands[4])->copy(),
                    (bands[5])->copy(),
                    (bands[6])->copy(),
                    (bands[7])->copy(),
                    (bands[8])->copy());

            return result;
        }
    };

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>::BandedMatrixQx(unsigned long size,
                    const DenseVector<DataType_> & ll,
                    const DenseVector<DataType_> & ld,
                    const DenseVector<DataType_> & lu,
                    const DenseVector<DataType_> & dl,
                    const DenseVector<DataType_> & dd,
                    const DenseVector<DataType_> & du,
                    const DenseVector<DataType_> & ul,
                    const DenseVector<DataType_> & ud,
                    const DenseVector<DataType_> & uu) :
        PrivateImplementationPattern<BandedMatrixQx<BandType_, DataType_>, Shared>(new Implementation<BandedMatrixQx<BandType_, DataType_> >(size))
    {
        CONTEXT("When creating BandedMatrixQx with initial bands:");
        ASSERT(size > 0, "size is zero!");


        if (ll.size() != size)
            throw VectorSizeDoesNotMatch(ll.size(), size);
        if (ld.size() != size)
            throw VectorSizeDoesNotMatch(ld.size(), size);
        if (lu.size() != size)
            throw VectorSizeDoesNotMatch(lu.size(), size);
        if (dl.size() != size)
            throw VectorSizeDoesNotMatch(dl.size(), size);
        if (dd.size() != size)
            throw VectorSizeDoesNotMatch(dd.size(), size);
        if (du.size() != size)
            throw VectorSizeDoesNotMatch(du.size(), size);
        if (ul.size() != size)
            throw VectorSizeDoesNotMatch(ul.size(), size);
        if (ud.size() != size)
            throw VectorSizeDoesNotMatch(ud.size(), size);
        if (uu.size() != size)
            throw VectorSizeDoesNotMatch(uu.size(), size);

        this->_imp->bands[LL].reset(new DenseVector<DataType_>(ll));
        this->_imp->bands[LD].reset(new DenseVector<DataType_>(ld));
        this->_imp->bands[LU].reset(new DenseVector<DataType_>(lu));
        this->_imp->bands[DL].reset(new DenseVector<DataType_>(dl));
        this->_imp->bands[DD].reset(new DenseVector<DataType_>(dd));
        this->_imp->bands[DU].reset(new DenseVector<DataType_>(du));
        this->_imp->bands[UL].reset(new DenseVector<DataType_>(ul));
        this->_imp->bands[UD].reset(new DenseVector<DataType_>(ud));
        this->_imp->bands[UU].reset(new DenseVector<DataType_>(uu));
    }

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>::BandedMatrixQx(BandedMatrix<DataType_> & src) :
        PrivateImplementationPattern<BandedMatrixQx<BandType_, DataType_>, Shared>(new Implementation<BandedMatrixQx<BandType_, DataType_> >(src))
    {
        CONTEXT("When creating BandedMatrixQx from BandedMatrix:");
    }

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>::BandedMatrixQx(SparseMatrixELL<DataType_> & src) :
        PrivateImplementationPattern<BandedMatrixQx<BandType_, DataType_>, Shared>(new Implementation<BandedMatrixQx<BandType_, DataType_> >(src))
    {
        CONTEXT("When creating BandedMatrixQx from SparseMatrixELL:");
    }

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>::BandedMatrixQx(const BandedMatrixQx<BandType_, DataType_> & other) :
        PrivateImplementationPattern<BandedMatrixQx<BandType_, DataType_>, Shared>(other._imp)
    {
    }

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>::~BandedMatrixQx()
    {
    }

    template <BandType BandType_, typename DataType_>
    unsigned long
    BandedMatrixQx<BandType_, DataType_>::columns() const
    {
        return this->_imp->size;
    }

    template <BandType BandType_, typename DataType_>
    unsigned long
    BandedMatrixQx<BandType_, DataType_>::rows() const
    {
        return this->_imp->size;
    }

    template <BandType BandType_, typename DataType_>
    unsigned long
    BandedMatrixQx<BandType_, DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <BandType BandType_, typename DataType_>
    const DataType_ &
    BandedMatrixQx<BandType_, DataType_>::operator() (unsigned long row, unsigned long column) const
    {
        CONTEXT("When retrieving BandedMatrixQx element, unassignable:");
        ASSERT(row < this->_imp->size, "row index '" + stringify(row) + "' is out of bounds!");
        ASSERT(column < this->_imp->size, "column index '" + stringify(column) + "' is out of bounds!");

        return (*this->_imp)(row, column);
    }

    template <BandType BandType_, typename DataType_>
    signed long
    BandedMatrixQx<BandType_, DataType_>::root() const
    {
        return this->_imp->root;
    }

    template <BandType BandType_, typename DataType_>
    unsigned long
    BandedMatrixQx<BandType_, DataType_>::bands() const
    {
        return this->_imp->band_count;
    }

    template <BandType BandType_, typename DataType_>
    DenseVector<DataType_> &
    BandedMatrixQx<BandType_, DataType_>::band(unsigned long index) const
    {
        CONTEXT("When retrieving band '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");

        return *this->_imp->bands[index];
    }

    template <BandType BandType_, typename DataType_>
    DenseVectorRange<DataType_>
    BandedMatrixQx<BandType_, DataType_>::band_range(unsigned long index) const
    {
        CONTEXT("When retrieving band range '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");

        return this->_imp->band_range(index);
    }

    template <BandType BandType_, typename DataType_>
    void BandedMatrixQx<BandType_, DataType_>::lock(LockMode mode) const
    {
        this->_imp->lock(mode);
    }

    template <BandType BandType_, typename DataType_>
            void BandedMatrixQx<BandType_, DataType_>::unlock(LockMode mode) const
    {
        this->_imp->unlock(mode);
    }

    template <BandType BandType_, typename DataType_>
    BandedMatrixQx<BandType_, DataType_>
    BandedMatrixQx<BandType_, DataType_>::copy() const
    {
        CONTEXT("When creating copy() of a BandedMatrixQx:");

        return this->_imp->copy();
    }

    template <BandType BandType_, typename DataType_>
    bool
    operator== (const BandedMatrixQx<BandType_, DataType_> & a, const BandedMatrixQx<BandType_, DataType_> & b)
    {
        if (a.size() != b.size())
            throw MatrixSizeDoesNotMatch(a.size(), b.size());

        bool result(true);

        for (unsigned long i(0) ; i < a.bands() ; ++i)
            result &= (a.band_range(i) == b.band_range(i));

        return result;
    }

    template <BandType BandType_, typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const BandedMatrixQx<BandType_, DataType_> & b)
    {
        lhs << "BandedMatrixQx" << std::endl << "{" << std::endl;
        for (unsigned long row(0) ; row < b.size() ; ++row)
        {
            for (unsigned long column(0) ; column < b.size() ; ++column)
            {
                lhs << " " << b(row, column);
            }
            lhs << std::endl;
        }
        lhs << "}" << std::endl;

        return lhs;
    }
}
#endif
