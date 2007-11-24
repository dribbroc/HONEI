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

#ifndef LIBLA_GUARD_ELEMENT_PRODUCT_HH
#define LIBLA_GUARD_ELEMENT_PRODUCT_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

#include <libla/dense_vector_range.hh>
#include <libutil/pool_task.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>

namespace honei
{
    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_ = tags::CPU> struct ElementProduct
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectorBases elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ;)
            {
                if (r.index() < l.index())
                {
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }
                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            return a;
            ///\todo: perhaps sparsify - in case l.index < r.index Write of 0 possible.
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVectorBase elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                *l *= b[l.index()];
            }
            return a;
            ///\todo: perhaps sparsify - if *b[l.index()] == 0 -> write of zero.
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix elements");

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
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and DenseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()); l != l_end ; ++l)
            {
                *l *= b[l.row()][l.column()];
            }

            return a; ///\todo: perhaps sparsify, dense_matrix[row][col] may be zero.
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                else if (r.index() < l.index())
                {
                    ++r;
                }

                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - in case l.index < r.index set to zero possible.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix elements");

            if (a.rows() != b.rows())
            {
                throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands());
            for (typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()),
                    l_end(a.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists()) //won't check l.exists() here, cause it is always created by Iterator.
                {
                    ++r;
                    continue;
                }

                ElementProduct<>::value(*l, *r);
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a DenseMatrix elements");

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
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a SparseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()); r != r_end ; )
            {
                while (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                *l *= *r;
                ++l; ++r;
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l_end(a.end_elements()) ;
                l != l_end ; ++l)
            {
                *l = DT1_(0);
            }
            /// \todo Left is complete dense at this point.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix and BandedMatrix elements");

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
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and BandedMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                while (r.index() < l.index())
                    ++r;

                *l *= *r;
                ++r; ++l;
            }

            return a;
        }

        /// \}

        #ifdef BENCHM
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            return result; 
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()); l != l_end ; ++l)
            {
                result.flops += 1;
                result.load += sizeof(DT1_) + sizeof(DT2_);
                result.store += 1;
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.size() * a.size();
            result.load = a.size() * a.size() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.size() * a.size() * sizeof(DT1_);
            return result;
        }
        #endif
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::Cell>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);
        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        /// \}
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::CPU::SSE>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & a, const DenseVector<double> & b);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };


    template <typename Tag_> struct MCElementProduct
    {
        template <typename DT1_, typename DT2_>
        static void value(SparseVector<DT1_> & a, const DenseVectorRange<DT2_> & b, unsigned long offset)
        {
            typename Vector<DT1_>::ElementIterator r(a.begin_non_zero_elements());
            r += offset;
            offset = r.index();
            unsigned long limit = r.index() + b.size();
            while (r.index() < limit)
            {
                *r *= b[r.index()-offset];
                ++r;
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of  DenseVector elements (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a,b);
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
                    TwoArgWrapper<ElementProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div, modulo+(i*div));
                    DenseVectorRange<DT2_> range_2(b, div, modulo+(i*div));
                    TwoArgWrapper<ElementProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
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
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVector elements (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long modulo = a.used_elements() % parts;
            unsigned long div = a.used_elements() / parts;
            if (div == 0) 
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                typename Vector<DT1_>::ElementIterator r(a.begin_non_zero_elements());
                unsigned long offset;
                for (int i(0); i < modulo; ++i) 
                {
                    offset = r.index();
                    r += div;
                    DenseVectorRange<DT2_> range(b, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCElementProduct<Tag_>, SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(a, range, (i*(div+1)));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT2_> range(b, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCElementProduct<Tag_>, SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(a, range, modulo + (i*div));
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
            CONTEXT("When calculating the product of BandedMatrix elements (MultiCore):");

            if (a.rows() != b.rows())
            {
                throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[2 * a.size()-1];
            int taskcount(0);
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands());
            for (typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()),
                    l_end(a.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists())
                {
                    ++r;
                    continue;
                }

                TwoArgWrapper< ElementProduct<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(*l, *r);
                pt[taskcount] = p->dispatch(mywrapper);
                ++r;
                ++taskcount;
            }
            for (int j = 0 ; j < taskcount ; ++j )
            {
                pt[j]->wait_on();
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix elements (MultiCore):");

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
                TwoArgWrapper< ElementProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and DenseMatrix elements (MultiCore):");

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
                TwoArgWrapper< ElementProduct<typename Tag_::DelegateTo>, SparseVector<DT1_>, const DenseVectorRange<DT2_> > mywrapper(a[i], b[i]);
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
            CONTEXT("When calculating the product of SparseMatrix elements (MulitCore):");

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
                TwoArgWrapper< ElementProduct<typename Tag_::DelegateTo>, SparseVector<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
            }
            return a;
        }
    };
    template <> struct ElementProduct <tags::CPU::MultiCore> : MCElementProduct <tags::CPU::MultiCore> {};
    template <> struct ElementProduct <tags::CPU::MultiCore::SSE> : MCElementProduct <tags::CPU::MultiCore::SSE> {};
}
#endif
