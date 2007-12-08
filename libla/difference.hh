/* vim: set sw=4 sts=4 et nofoldenable : */

/*
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

#ifndef LIBLA_GUARD_DIFFERENCE_HH
#define LIBLA_GUARD_DIFFERENCE_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/product.hh>
#include <libla/scale.hh>
#include <libla/sparse_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/vector.hh>
#include <libutil/tags.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Difference;

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Difference<tags::CPU>
    {
        /**
         * \name Matrix differences that return minuend
         * \{
         *
         * \brief Returns the the difference of two given matrices
         *
         * \param a The entity that is the minuend of the operation.
         * \param b The entity that is the subtrahend of the operation.
         *
         * \retval r Will modify the minuend a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a RowAccessMatrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                     r_end(b.end_non_zero_elements());
            for ( ; r != r_end ; ++r )
            {
                a[r.row()][r.column()] -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from SparseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
               if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = -(*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l -= *r;
                    ++l; ++r;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting BandedMatrix from BandedMatrix:");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                    Difference<>::value(*l, *r);
                else
                {
                    DenseVector<DT2_> band(r->copy());
                    Scale<>::value(DT1_(-1), band);
                    a.band(r.index() - a.size() + 1) = band;
                }
            }

            return a;
        }

        /// \}

        /**
         * \name Matrix differences that return subtrahend
         * \{
         *
         * \brief Returns the the difference of two given matrices.
         *
         * \param a The matrix that is the minuend of the operation.
         * \param b The matrix that is the subtrahend of the operation.
         *
         * \retval b Will modify the subtrahend b and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix:");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            Scale<tags::CPU>::value(-1,b);
            
            int middle_index(a.rows() -1);
            // If we are below the diagonal band, we start at Element index and go on until the last element.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (!vi.exists())
                    continue;
                unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                unsigned long i(0);
                for(typename Vector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                        c_end(vi->end_elements()) ; c != c_end ; ++c)
                {
                    b[start][i] += *c;
                    ++start, ++i;
                }
            }

            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
                    vi != vi_end ; ++vi)
            {
                if (!vi.exists())
                    continue;

                //Calculation of the element-index to stop in iteration!
                unsigned long offset(vi.index() - middle_index);
                unsigned long end(vi->size() - offset);
                unsigned long i(0);
                for(typename Vector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                        c_end(vi->element_at(end)) ; c != c_end ; ++c)
                {
                    b[i][offset] +=  *c;
                    ++offset, ++i;
                }
            }
            
            return b;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(const BandedMatrix<DT2_> & a, SparseMatrix<DT1_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix:");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }
            
            Scale<tags::CPU>::value(-1,b);
            
            int middle_index(a.rows() -1);
            // If we are below the diagonal band, we start at Element index and go on until the last element.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (!vi.exists())
                    continue;
                unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                unsigned long i(0);
                for(typename Vector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                        c_end(vi->end_elements()) ; c != c_end ; ++c)
                {
                    b[start][i] += *c;
                    ++start, ++i;
                }
            }

            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
                    vi != vi_end ; ++vi)
            {
                if (!vi.exists())
                    continue;

                //Calculation of the element-index to stop in iteration!
                unsigned long offset(vi.index() - middle_index);
                unsigned long end(vi->size() - offset);
                unsigned long i(0);
                for(typename Vector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                        c_end(vi->element_at(end)) ; c != c_end ; ++c)
                {
                    b[i][offset] +=  *c;
                    ++offset, ++i;
                }
            }

            return b;
        }

        /// \}

        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When subtracting DenseVectorBase from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l -= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.index()] = -(*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l -= *r;
                    ++l; ++r;
                }
            }
            return a;
            ///\todo: perhaps sparsify - i.e. substraction of 7 and 7 possible.
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                a[r.index()] -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT1_>::ConstElementIterator l(a.begin_elements());
            for (typename Vector<DT2_>::ElementIterator r(b.begin_elements()),
                    r_end(b.end_elements()) ; r != r_end ; ++r , ++l)
            {
                *r = *l - *r;
            }

            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT2_> & value(const SparseVector<DT1_> & a, DenseVector<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorRange<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorSlice<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }
        /// \}

        #ifdef BENCHM
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, SparseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            return result; 
        }
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            return result; 
        }

        #endif
    };

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::CPU::SSE>
    {
        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        static SparseMatrix<float> value(const BandedMatrix<float> &a, SparseMatrix<float> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static SparseMatrix<double> value(const BandedMatrix<double> &a, SparseMatrix<double> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static DenseMatrix<float> value(const BandedMatrix<float> &a, DenseMatrix<float> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static DenseMatrix<double> value(const BandedMatrix<double> &a, DenseMatrix<double> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        /// \}
    };

   /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::Cell>
    {
        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);
        static DenseVector<float> & value(const SparseVector<float> & a, DenseVector<float> & b);

        /// \}
    }; 

  //ANFANG

    template <typename Tag_> struct MCDifference
    { 
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & a, const SparseVector<DT2_> & b, unsigned long offset)
        {
            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            r += offset;
            offset = r.index();
            unsigned long limit = r.index() + a.size();
            while (r.index() < limit)
            {
                a[r.index()-offset] -= *r;
                ++r;
            }
        }
       
        template <typename DT1_, typename DT2_>
        static void value(const DenseVectorRange<DT1_> & a, SparseMatrix<DT2_> & b, unsigned long offset_y, unsigned long offset_x)
        {
            for(typename Vector<DT1_>::ConstElementIterator c(a.begin_elements()),
                        c_end(a.end_elements()) ; c != c_end ; ++c)
                {
                    b[offset_y][offset_x] += *c;
                    ++offset_y, ++offset_x;
                }

        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When substracting DenseVector from DenseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = a.size() % parts;
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                for (int i(0); i < modulo; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div+1, i*(div+1));
                    DenseVectorRange<DT2_> range_2(b, div+1, i*(div+1));
                    TwoArgWrapper<Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div, modulo+(i*div));
                    DenseVectorRange<DT2_> range_2(b, div, modulo+(i*div));
                    TwoArgWrapper<Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(0); i < parts; ++i)
                {
                    pt[i]->wait_on();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When substracting DenseVector from SparseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long modulo = b.used_elements() % parts;
            unsigned long div = b.used_elements() / parts;
            if (div == 0) 
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
                unsigned long offset;
                for (int i(0); i < modulo; ++i) 
                {
                    offset = r.index();
                    r += div;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, (i*(div+1)));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, modulo + (i*div));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i = 0; i < parts;  ++i)
                {
                    pt[i]->wait_on();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When substracting BandedMatrix from BandedMatrix (MultiCore):");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[2*a.rows()-1];
            int taskcount(0);
            typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                {
                    TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(*l, *r);
                    pt[taskcount] = p->dispatch(mywrapper);
                    ++taskcount; 
                }
                else
                {
                   DenseVector<DT2_> band(r->copy());
                   Scale<typename Tag_::DelegateTo>::value(DT1_(-1), band);
                   a.band(r.index()-a.size()+1) = band;
                }
            }
            for (unsigned long j = 0; j < taskcount; ++j)           
            {
                pt[j]->wait_on();
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        { 
            CONTEXT("When substracting DenseMatrix from DenseMatrix (MultiCore):");
            
            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        { 
            CONTEXT("When substracting DenseMatrix from SparseMatrix (MultiCore):");
            
            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        { 
            CONTEXT("When substracting SparseMatrix from SparseMatrix (MultiCore):");
            
            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, SparseVector<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(const BandedMatrix<DT2_> & a, SparseMatrix<DT1_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (MultiCore):");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }
            
            Scale<Tag_>::value(-1,b);

            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                return Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = a.size() % parts;
                ThreadPool * tp(ThreadPool::get_instance());
                std::list< PoolTask* > dispatched_tasks;
                Mutex mutex[parts];
                int middle_index(a.rows() -1);
                //if we are below the diagonal band
                for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(a.begin_bands()),
                     vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    unsigned long i(parts), offset(a.size());
                    unsigned long start(middle_index - vi.index());
                    while(i > modulo && offset-div > start)
                    {
                        --i;
                        offset-=div;
                        DenseVectorRange<DT2_> range_1(*vi, div, offset);
                        FourArgWrapper<MCDifference<Tag_>, const DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, (vi.index()+offset)-a.size());
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]); 
                        dispatched_tasks.push_back(tp->dispatch(func));
                    }
                    if (i == modulo)
                    {
                        while(i > 0 && offset-div-1 > start)
                        {
                            --i;
                            offset-=(div+1);
                            DenseVectorRange<DT2_> range_1(*vi, div+1, offset);
                            FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, (vi.index()+offset)-a.size());
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]); 
                            dispatched_tasks.push_back(tp->dispatch(func));
                        }
                    }
                    if (offset > start)
                    {
                        DenseVectorRange<DT2_> range_1(*vi, offset-start, start);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, start, 0);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i-1]); 
                        dispatched_tasks.push_back(tp->dispatch(func));
                    }
                }
                // If we are above or on the diagonal band
                for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(a.band_at(middle_index)), 
                     vi_end(a.end_bands()); vi != vi_end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    unsigned long i(0), offset(0);
                    unsigned long index(vi.index() - middle_index);
                    unsigned long end(vi->size() - index);
                    while ((i < modulo) && (offset+div+1 < end))
                    {
                        DenseVectorRange<DT2_> range_1(*vi, div+1, offset);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index + offset);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]); 
                        dispatched_tasks.push_back(tp->dispatch(func));
                        ++i;
                        offset+=div+1;
                    }
                    if (i == modulo)
                    {
                        while ((i < parts) && (offset+div  < end))
                        {
                            DenseVectorRange<DT2_> range_1(*vi, div, offset);
                            FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index+offset);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]); 
                            dispatched_tasks.push_back(tp->dispatch(func));
                            ++i;
                            offset+=div;
                        }
                    }
                    if (offset < end)
                    {
                        DenseVectorRange<DT2_> range_1(*vi, end - offset, offset);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index+offset);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]); 
                        dispatched_tasks.push_back(tp->dispatch(func));
                    }
                }
                
                while(! dispatched_tasks.empty())
                {
                    PoolTask * pt = dispatched_tasks.front();
                    dispatched_tasks.pop_front();
                    pt->wait_on();
                }
                return b;
            }
        }
    };
    template <> struct Difference <tags::CPU::MultiCore> : MCDifference <tags::CPU::MultiCore> {};
    template <> struct Difference <tags::CPU::MultiCore::SSE> : MCDifference <tags::CPU::MultiCore::SSE> {};

}
#endif
