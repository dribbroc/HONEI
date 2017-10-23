/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#pragma once
#ifndef LIBLA_GUARD_ALGORITHM_HH
#define LIBLA_GUARD_ALGORITHM_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/banded_matrix_qx.hh>
#include <honei/util/tags.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>
#include <honei/la/vector_error.hh>
#ifdef HONEI_GMP
#include <gmpxx.h>
#endif

#include <limits>

namespace honei
{
    /**
     * \{
     *
     * Copies data from a source container to a destination container.
     * The destination container data will reside in Tag_'s memory
     *
     * \ingroup grpalgorithm
     */

    template <typename Tag_, typename DT_> void copy(const DenseVectorContinuousBase<DT_> & source, DenseVectorContinuousBase<DT_> & dest)
    {
        CONTEXT("When copying elements from DenseVectorContinuousBase to DenseVectorContinuousBase:");
        ASSERT(source.elements() != dest.elements(),
                "trying to copy data from a DenseVectorContinuousBase to the very same DenseVectorContinuousBase!");

        if (source.size() != dest.size())
            throw VectorSizeDoesNotMatch(dest.size(), source.size());

        if (Tag_::tag_value == tags::tv_gpu_multi_core)
        {
            // TODO copy vektor parts direkt to the corresponding gpu memory
            /*DenseVectorRange<DT_> s1(source.range(source.size()/2, 0));
            DenseVectorRange<DT_> d1(dest.range(dest.size()/2, 0));
            DenseVectorRange<DT_> s2(source.range(source.size()/2 + source.size()%2, source.size()/2));
            DenseVectorRange<DT_> d2(dest.range(dest.size()/2 + dest.size()%2, dest.size()/2));
            cudaSumDVdouble task1(a1, b1, blocksize);
            cuda::GPUPool::instance()->enqueue(task2, 0)->wait();
            cuda::GPUPool::instance()->enqueue(task2, 1)->wait();*/
            dest.lock(lm_write_only);
            source.lock(lm_read_only);
            dest.unlock(lm_write_only);
            source.unlock(lm_read_only);
            TypeTraits<DT_>::copy(source.elements(), dest.elements(), dest.size());
        }
#ifdef HONEI_GMP
        else if (typeid(DT_) == typeid(mpf_class))
        {
            TypeTraits<DT_>::copy(source.elements(), dest.elements(), source.size());
        }
#endif
        else
        {
            MemoryArbiter::instance()->copy(Tag_::memory_value, source.memid(), source.address(),
                    dest.memid(), dest.address(), source.size() * sizeof(DT_));
        }
    }

    template <typename Tag_, typename DT_> void copy(const DenseVectorMPI<DT_> & source, DenseVectorMPI<DT_> & dest)
    {
        CONTEXT("When copying elements from DenseVectorMPI to DenseVectorMPI:");
        ASSERT(source.elements() != dest.elements(),
                "trying to copy data from a DenseVectorMPI to the very same DenseVectorMPI!");

        if (source.size() != dest.size())
            throw VectorSizeDoesNotMatch(dest.size(), source.size());

        /*if (Tag_::tag_value == tags::tv_gpu_multi_core)
        {
            dest.lock(lm_write_only);
            source.lock(lm_read_only);
            dest.unlock(lm_write_only);
            source.unlock(lm_read_only);
            TypeTraits<DT_>::copy(source.elements(), dest.elements(), dest.size());
        }
        else*/
        {
            MemoryArbiter::instance()->copy(Tag_::memory_value, source.vector().memid(), source.vector().address(),
                    dest.vector().memid(), dest.vector().address(), source.vector().size() * sizeof(DT_));
        }
    }

    template <typename IT_, typename DT_> void copy(const IT_ & begin, const IT_ & end,
            DenseVector<DT_> & dest)
    {
        CONTEXT("When copying elements from iterator range to DenseVector:");

        if (end.index() - begin.index() != dest.size())
            throw VectorSizeDoesNotMatch(dest.size(), end.index() - begin.index());

        typename DenseVector<DT_>::ConstElementIterator i(begin), i_end(end);
        for (typename DenseVector<DT_>::ElementIterator d(dest.begin_elements()) ; i != i_end ; ++i, ++d)
        {
            *d = *i;
        }
    }

    template <typename Tag_, typename DT_> void copy(const BandedMatrixQx<Q1Type, DT_> & source, BandedMatrixQx<Q1Type, DT_> & dest)
    {
        CONTEXT("When copying elements from BandedMatrixQ1 to BandedMatrixQ1:");

        if (source.size() != dest.size())
            throw VectorSizeDoesNotMatch(dest.size(), source.size());

        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(LL).memid(), source.band(LL).address(),
                dest.band(LL).memid(), dest.band(LL).address(), source.band(LL).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(LD).memid(), source.band(LD).address(),
                dest.band(LD).memid(), dest.band(LD).address(), source.band(LD).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(LU).memid(), source.band(LU).address(),
                dest.band(LU).memid(), dest.band(LU).address(), source.band(LU).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(DL).memid(), source.band(DL).address(),
                dest.band(DL).memid(), dest.band(DL).address(), source.band(DL).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(DD).memid(), source.band(DD).address(),
                dest.band(DD).memid(), dest.band(DD).address(), source.band(DD).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(DU).memid(), source.band(DU).address(),
                dest.band(DU).memid(), dest.band(DU).address(), source.band(DU).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(UL).memid(), source.band(UL).address(),
                dest.band(UL).memid(), dest.band(UL).address(), source.band(UL).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(UD).memid(), source.band(UD).address(),
                dest.band(UD).memid(), dest.band(UD).address(), source.band(UD).size() * sizeof(DT_));
        MemoryArbiter::instance()->copy(Tag_::memory_value, source.band(UU).memid(), source.band(UU).address(),
                dest.band(UU).memid(), dest.band(UU).address(), source.band(UU).size() * sizeof(DT_));
    }

    /// \}

    /**
     * \{
     *
     * Converts data from a source container to a destination container.
     *
     * \ingroup grpalgorithm
     */
    template <typename Tag_, typename DataType_> void convert(DenseVector<DataType_> & copy,
            const DenseVector<DataType_> & orig)
    {
        CONTEXT("When converting DenseVector to DenseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        MemoryArbiter::instance()->copy(Tag_::memory_value, orig.memid(), orig.address(),
                copy.memid(), copy.address(), orig.size() * sizeof(DataType_));
    }

    template <typename Tag_> void convert(DenseVector<double> & copy,
            const DenseVector<float> & orig)
    {
        CONTEXT("When converting DenseVector to DenseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        MemoryArbiter::instance()->convert_float_double(Tag_::memory_value, orig.memid(), orig.address(),
                copy.memid(), copy.address(), orig.size() * sizeof(float));
    }

    template <typename Tag_> void convert(DenseVector<float> & copy,
            const DenseVector<double> & orig)
    {
        CONTEXT("When converting DenseVector to DenseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        MemoryArbiter::instance()->convert_double_float(Tag_::memory_value, orig.memid(), orig.address(),
                copy.memid(), copy.address(), orig.size() * sizeof(double));
    }

#ifdef HONEI_GMP
    template <typename Tag_> void convert(DenseVector<mpf_class> & copy,
            const DenseVector<double> & orig)
    {
        CONTEXT("When converting DenseVector to DenseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        for (unsigned long i(0) ; i < copy.size() ; ++i)
        {
            mpf_class t(orig[i]);
            copy[i] = t;
        }
    }
#endif

    template <typename DataType_> void convert(DenseVectorBase<DataType_> & copy,
            const DenseVectorBase<DataType_> & orig)
    {
        CONTEXT("When converting DenseVectorBase to DenseVectorBase:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        copy = orig.copy();
    }

    template <typename DataType_, typename OrigType_> void convert(DenseVectorBase<DataType_> & copy,
            const DenseVectorBase<OrigType_> & orig)
    {
        CONTEXT("When converting DenseVectorBase to DenseVectorBase:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        copy.lock(lm_write_only);
        orig.lock(lm_read_only);
        typename DenseVector<DataType_>::ElementIterator f(copy.begin_elements());
        for (typename DenseVector<OrigType_>::ConstElementIterator i(orig.begin_elements()),
                i_end(orig.end_elements()) ; i != i_end ; ++i, ++f)
        {
            *f = *i;
        }
        copy.unlock(lm_write_only);
        orig.unlock(lm_read_only);
    }

    template <typename DataType_> void convert(SparseVector<DataType_> & copy,
            const SparseVector<DataType_> & orig)
    {
        CONTEXT("When converting SparseVector to SparseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        copy = orig.copy();
    }

    template <typename DataType_, typename OrigType_> void convert(SparseVector<DataType_> & copy,
            const SparseVector<OrigType_> & orig)
    {
        CONTEXT("When converting SparseVector to SparseVector:");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        for (typename SparseVector<OrigType_>::NonZeroConstElementIterator i(orig.begin_non_zero_elements()),
                i_end(orig.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            copy[i.index()] = *i;
        }
    }

    template <typename DataType_> void convert(DenseMatrix<DataType_> & copy,
            const DenseMatrix<DataType_> & orig)
    {
        CONTEXT("When converting DenseMatrix to DenseMatrix:");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        copy = orig.copy();
    }

    template <typename DataType_, typename OrigType_> void convert(DenseMatrix<DataType_> & copy,
            const DenseMatrix<OrigType_> & orig)
    {
        CONTEXT("When converting DenseMatrix to DenseMatrix:");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        copy.lock(lm_write_only);
        orig.lock(lm_read_only);
        TypeTraits<OrigType_>::convert(copy.elements(), orig.elements(), orig.columns() * orig.rows());
        copy.unlock(lm_write_only);
        orig.unlock(lm_read_only);
    }

    template <typename DataType_> void convert(SparseMatrix<DataType_> & copy,
            const SparseMatrix<DataType_> & orig)
    {
        CONTEXT("When converting SparseMatrix to SparseMatrix:");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        copy = orig.copy();
    }

    template <typename DataType_, typename OrigType_> void convert(SparseMatrix<DataType_> & copy,
            const SparseMatrix<OrigType_> & orig)
    {
        CONTEXT("When converting SparseMatrix to SparseMatrix:");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        for (unsigned long i(0) ; i < orig.rows() ; ++i)
        {
            convert(copy[i], orig[i]);
        }
    }

    template <typename DataType_> void convert(BandedMatrix<DataType_> & copy,
            const BandedMatrix<DataType_> & orig)
    {
        CONTEXT("When converting BandedMatrix to BandedMatrix:");

        if (copy.size() != orig.size())
            throw MatrixSizeDoesNotMatch(orig.size(), copy.size());

        copy = orig.copy();
    }


    template <typename DataType_, typename OrigType_> void convert(BandedMatrix<DataType_> & copy,
            const BandedMatrix<OrigType_> & orig)
    {
        CONTEXT("When converting BandedMatrix to BandedMatrix:");

        if (copy.size() != orig.size())
            throw MatrixSizeDoesNotMatch(orig.size(), copy.size());

        for (typename BandedMatrix<OrigType_>::ConstBandIterator band(orig.begin_non_zero_bands()), band_end(orig.end_non_zero_bands()) ;
                band != band_end ; ++band)
        {
            convert(copy.band_unsigned(band.index()), *band);
        }
    }

    /// \}

    /**
     * \{
     *
     * Fills a container with a given prototype.
     *
     * \ingroup grpalgorithm
     */

    template <typename Tag_>
    void fill(DenseVectorContinuousBase<float> & dest, const float & proto = float(0))
    {
        CONTEXT("When filling DenseVectorContinuousBase with '" + stringify(proto) + "':");

        //TypeTraits<DT_>::fill(dest.elements(), dest.size(), proto);
        MemoryArbiter::instance()->fill(Tag_::memory_value, dest.memid(), dest.address(),
                dest.size() * sizeof(float), proto);
    }

    template <typename Tag_>
    void fill(DenseVectorContinuousBase<double> & dest, const double & proto = double(0))
    {
        CONTEXT("When filling DenseVectorContinuousBase with '" + stringify(proto) + "':");

        //TypeTraits<DT_>::fill(dest.elements(), dest.size(), proto);
        MemoryArbiter::instance()->fill(Tag_::memory_value, dest.memid(), dest.address(),
                dest.size() * sizeof(double), proto);
    }

    template <typename Tag_, typename DT_>
    void fill(DenseVectorMPI<DT_> & dest, const DT_ & proto = DT_(0))
    {
        CONTEXT("When filling DenseVectorMPI with '" + stringify(proto) + "':");

        MemoryArbiter::instance()->fill(Tag_::memory_value, dest.vector().memid(), dest.vector().address(),
                dest.vector().size() * sizeof(DT_), proto);
    }

    template <typename Tag_>
    void fill(const DenseMatrix<float> & dest, const float & proto = float(0))
    {
        CONTEXT("When filling DenseMatrix with '" + stringify(proto) + "':");

        //TypeTraits<DT_>::fill(dest.elements(), dest.rows() * dest.columns(), proto);
        MemoryArbiter::instance()->fill(Tag_::memory_value, dest.memid(), dest.address(),
                dest.size() * sizeof(float), proto);
    }

    template <typename Tag_>
    void fill(const DenseMatrix<double> & dest, const double & proto = double(0))
    {
        CONTEXT("When filling DenseMatrix with '" + stringify(proto) + "':");

        //TypeTraits<DT_>::fill(dest.elements(), dest.rows() * dest.columns(), proto);
        MemoryArbiter::instance()->fill(Tag_::memory_value, dest.memid(), dest.address(),
                dest.size() * sizeof(double), proto);
    }

    template <typename Tag, typename DT_> void fill(const DenseVectorContinuousBase<DT_> & dest, const DT_ & proto = DT_(0))
    {
        CONTEXT("When filling DenseVectorContinuousBase with '" + stringify(proto) + "':");

        dest.lock(lm_write_only);
        TypeTraits<DT_>::fill(dest.elements(), dest.size(), proto);
        dest.unlock(lm_write_only);
    }

    template <typename Tag, typename DT_> void fill(const DenseMatrix<DT_> & dest, const DT_ & proto = DT_(0))
    {
        CONTEXT("When filling DenseMatrix with '" + stringify(proto) + "':");

        dest.lock(lm_write_only);
        TypeTraits<DT_>::fill(dest.elements(), dest.rows() * dest.columns(), proto);
        dest.unlock(lm_write_only);
    }

    template <typename IT_, typename DT_> void fill(const IT_ & begin, const IT_ & end, const DT_ & proto = DT_(0))
    {
        CONTEXT("When filling elements of iterator range with '" + stringify(proto) + "':");

        for (IT_ i(begin), i_end(end) ; i != i_end ; ++i)
        {
            *i = proto;
        }
    }

    /// \}

}

#endif
