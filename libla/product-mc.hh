/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 David Gies <david-gies@gmx.de>
 * Copyright (c) 2007 Joachim Messer <joachim.messer@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_PRODUCT_MC_HH
#define LIBLA_GUARD_PRODUCT_MC_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/dot_product.hh>
#include <libla/element_product.hh>
#include <libla/matrix_error.hh>
#include <libla/scaled_sum.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/sum.hh>
#include <libutil/pool_task.hh>
#include <libutil/tags.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>

#include <tr1/functional>
#include <cmath>
#include <iostream>

namespace honei
{
    template <typename Tag_> struct Product;

    /**
     * \brief Product of two entities.
     *
     * MatrixProduct is the class template for the product operation
     * \f[
     *     \texttt{Product}(a, b): \quad c \leftarrow a * b,
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_> struct MCProduct
    {
        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseVector(Base) (MultiCore):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows(), DT1_(0));
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];

            for (unsigned long i = 0; i < a.rows(); ++i)
            {

                ResultTwoArgWrapper< DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>,
                    const DenseVector<DT2_> > mywrapper(*l, a[i], b);
                pt[i] = p->dispatch(mywrapper);

                ++l;
            }

            for (unsigned long i = 0; i < a.rows();  ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix (MultiCore):");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask   * pt[a.rows() * b.columns()];


            DenseMatrix<DT1_> result(a.rows(), b.columns());

            for (unsigned int s(0) ; s < a.rows() ; ++s)
            {
                const DenseVectorRange<DT1_> a_row(a[s]);
                for (unsigned int t(0); t < b.columns() ; ++t)
                {
                    const DenseVectorSlice<DT2_> b_column(b.column(t));
                    //result[s][t] = DotProduct<>::value(b_column, a_row);
                    ResultTwoArgWrapper< DotProduct<>, DT1_, const DenseVectorRange<DT1_>,
                        const DenseVectorSlice<DT2_> > mywrapper(result[s][t], a_row, b_column);
                    pt[s * b.columns() + t] = p->dispatch(mywrapper);
                }

            }

            for (unsigned long i = 0; (i < a.rows() * b.columns());  ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return result;
        }


       // HelpFunction for BandedMatrix * DenseVector MultiCore
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & result, const DenseVectorRange<DT1_> & a, const DenseVectorRange<DT2_> & b, Mutex * mutex)
        {
            DenseVector<DT1_> temp_result(a.copy());
            ElementProduct<typename Tag_::DelegateTo>::value(temp_result, b);
            Lock l(*mutex);
            Sum<typename Tag_::DelegateTo>::value(result, temp_result);
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const BandedMatrix<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with DenseVector(MultiCore):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows(), DT1_(0));
            unsigned long parts(8);
            unsigned long div = b.size() / parts;
            if (div == 0)
            {
                return Product<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = b.size() % parts;
                ThreadPool * tp(ThreadPool::get_instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                Mutex mutex[parts];
                int middle_index(a.rows() -1);

                for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                        vi != vi_end ; ++vi)
                {
                    // If we are above or on the diagonal band
                    if (middle_index < vi.index())
                    {
                        if (!vi.exists())
                            continue;
                        unsigned long i(0), offset(0);
                        unsigned long end(vi->size() - (vi.index() - middle_index));
                        while ((i < modulo) && (offset+div+1 < end))
                        {
                            DenseVectorRange<DT1_> range_1(*vi, div+1, offset);
                            DenseVectorRange<DT2_> range_2(b, div+1, offset + vi.index() - middle_index);
                            DenseVectorRange<DT1_> res_range(result, div+1, offset);
                            ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                            ++i;
                            offset+=div+1;
                        }
                        if (i == modulo)
                        {
                            while ((i < parts) && (offset+div  < end))
                            {
                                DenseVectorRange<DT1_> range_1(*vi, div, offset);
                                DenseVectorRange<DT2_> range_2(b, div, offset + vi.index() - middle_index);
                                DenseVectorRange<DT1_> res_range(result, div, offset);
                                ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                    DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                                std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                                std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                                dispatched_tasks.push_back(ptr);
                                ++i;
                                offset+=div;
                            }
                        }
                        if (offset < end)
                        {
                            DenseVectorRange<DT1_> range_1(*vi, end - offset, offset);
                            DenseVectorRange<DT2_> range_2(b, end - offset, offset + vi.index()-middle_index);
                            DenseVectorRange<DT1_> res_range(result, end - offset, offset);
                            ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                        }
                    }
                    //if we are below the diagonal band
                    else
                    {
                        if (!vi.exists())
                            continue;
                        unsigned long i(parts), offset(b.size());
                        unsigned long start(middle_index - vi.index());
                        while(i > modulo && offset-div > start)
                        {
                            --i;
                            offset-=div;
                            DenseVectorRange<DT1_> range_1(*vi, div, offset);
                            DenseVectorRange<DT2_> range_2(b, div, offset - (middle_index - vi.index()));
                            DenseVectorRange<DT1_> res_range(result, div, offset);
                            ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                        }
                        if (i == modulo)
                        {
                            while(i > 0 && offset-div-1 > start)
                            {
                                --i;
                                offset-=(div+1);
                                DenseVectorRange<DT1_> range_1(*vi, div+1, offset);
                                DenseVectorRange<DT2_> range_2(b, div+1, offset - (middle_index - vi.index()));
                                DenseVectorRange<DT1_> res_range(result, div+1, offset);
                                ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                    DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                                std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                                std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                                dispatched_tasks.push_back(ptr);
                            }
                        }
                        if (offset > start)
                        {
                            DenseVectorRange<DT1_> range_1(*vi, offset-start, start);
                            DenseVectorRange<DT2_> range_2(b, offset-start, 0);
                            DenseVectorRange<DT1_> res_range(result, offset-start, start);
                            ThreeArgWrapper<ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                                DenseVectorRange<DT1_>, DenseVectorRange<DT2_> > mywrapper(res_range, range_1, range_2);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i-1]);
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                        }
                    }
                }

                while(! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
                return result;
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseMatrix (MultiCore):");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            Mutex mutex[a.rows()];

            for (typename SparseMatrix<DT1_>::ConstRowIterator i(a.begin_non_zero_rows()), i_end(a.end_non_zero_rows()) ;
                    i != i_end ; ++i)
            {
                for (typename Vector<DT1_>::ConstElementIterator vi((*i).begin_non_zero_elements()), vi_end((*i).end_non_zero_elements()) ;
                        vi != vi_end ; ++vi)
                {
                    ThreeArgWrapper< ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_>, const DT1_ >
                        wrapper(result[i.index()], b[vi.index()], *vi);
                    std::tr1::function<void ()> func = std::tr1::bind(wrapper, &mutex[i.index()]);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                    dispatched_tasks.push_back(ptr);
                }
            }

            while (! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return result;
        }
    };
}
#endif
