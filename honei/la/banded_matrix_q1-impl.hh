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

#ifndef LIBLA_GUARD_BANDED_MATRIX_Q1_IMPL_HH
#define LIBLA_GUARD_BANDED_MATRIX_Q1_IMPL_HH 1

#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/assertion.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

namespace honei
{
    template <typename DataType_> struct Implementation<BandedMatrixQ1<DataType_> >
    {
        /// Array of pointers to our band-data.
        SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > bands;

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
            size(size),
            root((unsigned long)sqrt(size)),
            zero_vector(size, DataType_(0))
        {
        }
    };

    template <typename DataType_>
    BandedMatrixQ1<DataType_>::BandedMatrixQ1(unsigned long size,
                    const DenseVector<DataType_> & ll,
                    const DenseVector<DataType_> & ld,
                    const DenseVector<DataType_> & lu,
                    const DenseVector<DataType_> & dl,
                    const DenseVector<DataType_> & dd,
                    const DenseVector<DataType_> & du,
                    const DenseVector<DataType_> & ul,
                    const DenseVector<DataType_> & ud,
                    const DenseVector<DataType_> & uu) :
        PrivateImplementationPattern<BandedMatrixQ1<DataType_>, Shared>(new Implementation<BandedMatrixQ1<DataType_> >(size))
    {
        CONTEXT("When creating BandedMatrixQ1 with initial bands:");
        ASSERT(size > 0, "size is zero!");

        unsigned long target;
        for (unsigned long level(0) ; level <= 50 ; ++level)
        {
            target = (unsigned long)pow((pow(2, level) + 1), 2);
            if (target >= size)
                break;
        }
        ASSERT(size == target, "size of " + stringify(size) + " matches not ((2^L) +1)^2, L integer.");

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

    template <typename DataType_>
    BandedMatrixQ1<DataType_>::BandedMatrixQ1(BandedMatrix<DataType_> & src) :
        PrivateImplementationPattern<BandedMatrixQ1<DataType_>, Shared>(new Implementation<BandedMatrixQ1<DataType_> >(src.size()))
    {
        CONTEXT("When creating BandedMatrixQ1 from BandedMatrix:");
        ASSERT(src.size() > 0, "size is zero!");

        unsigned long target;
        for (unsigned long level(1) ; level <= 50 ; ++level)
        {
            target = (unsigned long)pow((pow(2, level) + 1), 2);
            if (target >= src.size())
                break;
        }
        ASSERT(src.size() == target, "size matches not ((2^L) +1)^2, L integer.");

        this->_imp->bands[LL].reset(new DenseVector<DataType_>(src.band(-this->_imp->root - 1)));
        this->_imp->bands[LD].reset(new DenseVector<DataType_>(src.band(-this->_imp->root)));
        this->_imp->bands[LU].reset(new DenseVector<DataType_>(src.band(-this->_imp->root + 1)));
        this->_imp->bands[DL].reset(new DenseVector<DataType_>(src.band(-1)));
        this->_imp->bands[DD].reset(new DenseVector<DataType_>(src.band(0)));
        this->_imp->bands[DU].reset(new DenseVector<DataType_>(src.band(1)));
        this->_imp->bands[UL].reset(new DenseVector<DataType_>(src.band(this->_imp->root - 1)));
        this->_imp->bands[UD].reset(new DenseVector<DataType_>(src.band(this->_imp->root)));
        this->_imp->bands[UU].reset(new DenseVector<DataType_>(src.band(this->_imp->root + 1)));
    }

    template <typename DataType_>
    BandedMatrixQ1<DataType_>::~BandedMatrixQ1()
    {
    }

    template <typename DataType_>
    unsigned long
    BandedMatrixQ1<DataType_>::columns() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long
    BandedMatrixQ1<DataType_>::rows() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long
    BandedMatrixQ1<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    const DataType_ &
    BandedMatrixQ1<DataType_>::operator() (unsigned long row, unsigned long column) const
    {
        CONTEXT("When retrieving BandedMatrixQ1 element, unassignable:");
        ASSERT(row < this->_imp->size && row >= 0, "row index '" + stringify(row) + "' is out of bounds!");
        ASSERT(column < this->_imp->size && column >= 0, "column index '" + stringify(column) + "' is out of bounds!");

        long index(column - row);
        long root(this->_imp->root);

        if (index == -root - 1)
            return (*this->_imp->bands[LL])[row];
        else if (index == -root)
            return (*this->_imp->bands[LD])[row];
        else if (index == -root + 1)
            return (*this->_imp->bands[LU])[row];
        else if (index == -1)
            return (*this->_imp->bands[DL])[row];
        else if (index == 0)
            return (*this->_imp->bands[DD])[row];
        else if (index == 1)
            return (*this->_imp->bands[DU])[row];
        else if (index == root - 1)
            return (*this->_imp->bands[UL])[row];
        else if (index == root)
            return (*this->_imp->bands[UD])[row];
        else if (index == root + 1)
            return (*this->_imp->bands[UU])[row];
        else
            return this->_imp->zero_element;
    }

    template <typename DataType_>
    signed long
    BandedMatrixQ1<DataType_>::root() const
    {
        return this->_imp->root;
    }

    template <typename DataType_>
    DenseVector<DataType_> &
    BandedMatrixQ1<DataType_>::band(Q1BandIndex index) const
    {
        CONTEXT("When retrieving band '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");

        return *this->_imp->bands[index];
    }

    template <typename DataType_>
    DenseVectorRange<DataType_>
    BandedMatrixQ1<DataType_>::band_range(Q1BandIndex index) const
    {
        CONTEXT("When retrieving band range '" + stringify(index) + "' of matrix of size '"
                + stringify(this->_imp->size) + "':");

        switch(index)
        {
            case LL:
                return DenseVectorRange<DataType_>(*this->_imp->bands[LL], this->_imp->size - this->_imp->root - 1, this->_imp->root + 1);
                break;
            case LD:
                return DenseVectorRange<DataType_>(*this->_imp->bands[LD], this->_imp->size - this->_imp->root, this->_imp->root);
                break;
            case LU:
                return DenseVectorRange<DataType_>(*this->_imp->bands[LU], this->_imp->size - this->_imp->root + 1, this->_imp->root - 1 );
                break;
            case DL:
                return DenseVectorRange<DataType_>(*this->_imp->bands[DL], this->_imp->size - 1, 1);
                break;
            case DD:
                return DenseVectorRange<DataType_>(*this->_imp->bands[DD], this->_imp->size, 0);
                break;
            case DU:
                return DenseVectorRange<DataType_>(*this->_imp->bands[DU], this->_imp->size - 1, 0);
                break;
            case UL:
                return DenseVectorRange<DataType_>(*this->_imp->bands[UL], this->_imp->size - this->_imp->root + 1, 0);
                break;
            case UD:
                return DenseVectorRange<DataType_>(*this->_imp->bands[UD], this->_imp->size - this->_imp->root, 0);
                break;
            case UU:
                return DenseVectorRange<DataType_>(*this->_imp->bands[UU], this->_imp->size - this->_imp->root - 1, 0);
                break;
            default:
                return DenseVectorRange<DataType_>(this->_imp->zero_vector, this->_imp->size, 0);
        }
    }

    template <typename DataType_>
    void * BandedMatrixQ1<DataType_>::lock(LockMode mode, tags::TagValue memory) const
    {
        this->_imp->bands[LL]->lock(mode, memory);
        this->_imp->bands[LD]->lock(mode, memory);
        this->_imp->bands[LU]->lock(mode, memory);
        this->_imp->bands[DL]->lock(mode, memory);
        this->_imp->bands[DD]->lock(mode, memory);
        this->_imp->bands[DU]->lock(mode, memory);
        this->_imp->bands[UL]->lock(mode, memory);
        this->_imp->bands[UD]->lock(mode, memory);
        this->_imp->bands[UU]->lock(mode, memory);
    }

    template <typename DataType_>
            void BandedMatrixQ1<DataType_>::unlock(LockMode mode) const
    {
        this->_imp->bands[LL]->unlock(mode);
        this->_imp->bands[LD]->unlock(mode);
        this->_imp->bands[LU]->unlock(mode);
        this->_imp->bands[DL]->unlock(mode);
        this->_imp->bands[DD]->unlock(mode);
        this->_imp->bands[DU]->unlock(mode);
        this->_imp->bands[UL]->unlock(mode);
        this->_imp->bands[UD]->unlock(mode);
        this->_imp->bands[UU]->unlock(mode);
    }

    template <typename DataType_>
    BandedMatrixQ1<DataType_>
    BandedMatrixQ1<DataType_>::copy() const
    {
        CONTEXT("When creating copy() of a BandedMatrixQ1:");
        BandedMatrixQ1 result(this->_imp->size,
                (this->_imp->bands[0])->copy(),
                (this->_imp->bands[1])->copy(),
                (this->_imp->bands[2])->copy(),
                (this->_imp->bands[3])->copy(),
                (this->_imp->bands[4])->copy(),
                (this->_imp->bands[5])->copy(),
                (this->_imp->bands[6])->copy(),
                (this->_imp->bands[7])->copy(),
                (this->_imp->bands[8])->copy());

        return result;
    }

    template <typename DataType_>
    bool
    operator== (const BandedMatrixQ1<DataType_> & a, const BandedMatrixQ1<DataType_> & b)
    {
        if (a.size() != b.size())
            throw MatrixSizeDoesNotMatch(a.size(), b.size());

        bool result(true);

        result &= (a.band_range(LL) == b.band_range(LL));
        result &= (a.band_range(LD) == b.band_range(LD));
        result &= (a.band_range(LU) == b.band_range(LU));
        result &= (a.band_range(DL) == b.band_range(DL));
        result &= (a.band_range(DD) == b.band_range(DD));
        result &= (a.band_range(DU) == b.band_range(DU));
        result &= (a.band_range(UL) == b.band_range(UL));
        result &= (a.band_range(UD) == b.band_range(UD));
        result &= (a.band_range(UU) == b.band_range(UU));

        return result;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const BandedMatrixQ1<DataType_> & b)
    {
        lhs << "BandedMatrixQ1" << std::endl << "{" << std::endl;
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
