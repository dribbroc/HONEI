/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_ELEMENT_INVERSE_HH
#define LIBLA_GUARD_ELEMENT_INVERSE_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/tags.hh>

#include <algorithm>

namespace honei
{
    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_ = tags::CPU> struct ElementInverse
    {
        /**
         * \name Element inversions
         * \{
         *
         * \brief Returns the inverse values of all of an entity's elements.
         *
         * \param x The entity whose elements' inverse values shall be computed.
         *
         * \retval x Will modify the entity x and return it.
         */

        template <typename DataType_>
        static DenseVectorBase<DataType_> & value(DenseVectorBase<DataType_> & x)
        {
            CONTEXT("When calculating the inverse DenseVectorBase elements:");

            for (typename DenseVectorBase<DataType_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                 if (*l == DataType_(0))
                    continue;

                *l = DataType_(1) / *l;
            }
            return x;
        }

        template <typename DataType_>
        static SparseVector<DataType_> & value(SparseVector<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseVector elements:");


            for (typename Vector<DataType_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                  if (*l == DataType_(0))
                    continue;

                  *l = DataType_(1) / *l;
            }
            return x;
        }

        template <typename DataType_>
        static DenseMatrix<DataType_> & value(DenseMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseMatrix elements:");

            for (typename DenseMatrix<DataType_>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i == DataType_(0))
                    continue;

                *i = DataType_(1) / *i;
            }

            return x;
        }

        template <typename DataType_>
        static SparseMatrix<DataType_> & value(SparseMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseMatrix elements:");

            for (typename MutableMatrix<DataType_>::ElementIterator i(x.begin_non_zero_elements()),
                    i_end(x.end_non_zero_elements()) ; i != i_end ; ++i)
            {
               if (*i == DataType_(0))
                  continue;

               *i = DataType_(1) / *i;
            }

            return x;
        }

        template <typename DataType_>
        static BandedMatrix<DataType_> & value(BandedMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of BandedVector elements:");

            for (typename BandedMatrix<DataType_>::BandIterator i(x.begin_non_zero_bands()),
                    i_end(x.end_non_zero_bands()) ; i != i_end ; ++i)
            {
                DenseVector<DataType_> band = *i;
                unsigned long middle_index(x.rows() - 1);
                // If we are above or on the diagonal band, we start at Element 0 and go on
                // until Element band_size-band_index.
                if (i.index() >= middle_index)
                {
                    //Calculation of the element-index to stop in iteration!
                    unsigned long end = band.size() - (i.index() - middle_index);

                    for (typename DenseVector<DataType_>::ElementIterator b(band.begin_elements()),
                            b_end(band.element_at(end)) ; b < b_end ; ++b)
                    {
                        if (*b == DataType_(0))
                            continue;

                        *b = DataType_(1) / *b;
                    }
                }
                else
                {
                    //Calculation of the element-index to start in iteration!
                    unsigned long start = middle_index - i.index();
                    for (typename DenseVector<DataType_>::ElementIterator b(band.element_at(start)),
                            b_end(band.end_elements()) ; b < b_end ; ++b)
                    {
                        if (*b == DataType_(0))
                            continue;

                        *b = DataType_(1) / *b;
                    }
                }
            }

            return x;
        }

        template <typename DT_>
        static inline DenseVector<DT_> & value(DenseVector<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            ElementInverse<>::value(temp);
            return x;
        }

        template <typename DT_>
        static inline DenseVectorRange<DT_> & value(DenseVectorRange<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            ElementInverse<>::value(temp);
            return x;
        }

        template <typename DT_>
        static inline DenseVectorSlice<DT_> & value(DenseVectorSlice<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            ElementInverse<>::value(temp);
            return x;
        }

        /// \}

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * sizeof(DT1_);
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            result.size.push_back(result.flops);
            return result; 
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseVector<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = a.size() * sizeof(DT1_);
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            return result; 
        }
    };

    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct ElementInverse<tags::GPU::CUDA>
    {
        /**
         * \name Element inversions
         * \{
         *
         * \brief Returns the inverse values of all of an entity's elements.
         *
         * \param x The entity whose elements' inverse values shall be computed.
         *
         * \retval x Will modify the entity x and return it.
         */
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x);

        static DenseMatrix<float> & value(DenseMatrix<float> & x);

        /// \}
    };

    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct ElementInverse<tags::CPU::SSE>
    {
        /**
         * \name Element inversions
         * \{
         *
         * \brief Returns the inverse values of all of an entity's elements.
         *
         * \param x The entity whose elements' inverse values shall be computed.
         *
         * \retval x Will modify the entity x and return it.
         */
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x);

        static DenseMatrix<float> & value(DenseMatrix<float> & x);

        static DenseMatrix<double> & value(DenseMatrix<double> & x);

        static SparseVector<float> & value(SparseVector<float> & x);

        static SparseVector<double> & value(SparseVector<double> & x);

        static SparseMatrix<float> & value(SparseMatrix<float> & x);

        static SparseMatrix<double> & value(SparseMatrix<double> & x);
        /// \}
    };

    namespace mc
    {
        template <typename Tag_> struct ElementInverse :
            public honei::ElementInverse<typename Tag_::DelegateTo>
        {
        };
    }

    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */

    template <> struct ElementInverse<tags::CPU::MultiCore> :
        public mc::ElementInverse<tags::CPU::MultiCore>
    {
    };

    template <> struct ElementInverse<tags::CPU::MultiCore::SSE> :
        public mc::ElementInverse<tags::CPU::MultiCore::SSE>
    {
    };

    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementInverse<tags::Cell>
    {
        /**
         * \name Element inversions
         * \{
         *
         * \brief Returns the inverse values of all of an entity's elements.
         *
         * \param a The entity whose elements' inverse values shall be computed.
         *
         * \retval a Will modify the entity x and return it.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a);
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a);

        static DenseMatrix<double> & value(DenseMatrix<double> & a)
        {
            CONTEXT("When forwarding ElementInverse tags::Cell to tags::CPU:");
            ElementInverse<tags::CPU>::value(a);
            return a;
        }
        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a)
        {
            CONTEXT("When forwarding ElementInverse tags::Cell to tags::CPU:");
            ElementInverse<tags::CPU>::value(a);
            return a;
        }
    };
}
#endif
