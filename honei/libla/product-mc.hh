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

#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/dot_product.hh>
#include <honei/libla/element_product.hh>
#include <honei/libla/matrix_error.hh>
#include <honei/libla/scaled_sum.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/sum.hh>
#include <honei/libla/scaled_sum.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/tags.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>

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
        // Help function for DenseMatrix * DenseVector MultiCore.
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & result, const DenseMatrix<DT1_> & a,
                const DenseVectorRange<DT2_> & b, const unsigned long offset, const unsigned long size)
        {
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());

            for (unsigned long i(offset) ; i < offset + size ; ++i)
            {
                *l = DotProduct<Tag_>::value(a[i], b);
                ++l;
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseVectorContinuousBase (MultiCore):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::product[DM,DVCB]::min-part-size", 256));
            unsigned long overall_size(a.rows() * a.columns());

            if (overall_size < (min_part_size << 1))
            {
                return Product<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                DenseVector<DT1_> result(a.rows(), DT1_(0));
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::product[DM,DVCB]::max-count", num_threads));
                if (min_part_size > a.columns())
                {
                    min_part_size = min_part_size / a.columns() + (min_part_size % a.columns())? 1 : 0;
                }
                else
                {
                    min_part_size = 1;
                }
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 1,  a.rows(), PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    DenseVectorRange<DT1_> res_range(result, part_size, offset);
                    FiveArgWrapper< MCProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                        const DenseMatrix<DT1_>, const DenseVectorRange<DT2_>, const unsigned long,
                        const unsigned long >
                        mywrapper(res_range, a, b.range(b.size(), 0), offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }

                return result;
            }
        }
/*
        // Help function for DenseMatrix * DenseMatrix MultiCore.
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & a, const DenseMatrix<DT2_> & b, const DenseVectorRange<DT1_> & c)
        {
            CONTEXT("When partial multiplying DenseMatrix with DenseMatrix:");
            
            for (unsigned long j(0) ; j < b.rows() ; ++j)
            {
                ScaledSum<Tag_>::value(a, b[j], DT2_(c[j]));
            }
        }

        
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix (MultiCore):");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            ThreadPool * p(ThreadPool::instance());
            PoolTask   * pt[a.rows()];

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));
            
            for (unsigned long i(0) ; i < a.rows() ; ++i)
            {
                ThreeArgWrapper< MCProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                     const DenseMatrix<DT2_>, const DenseVectorRange<DT1_> > mywrapper(result[i], b, a[i]);
                pt[i] = p->dispatch(mywrapper);
            }

            for (unsigned long i = 0; (i < a.rows());  ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            } 
            
            return result;
        }

*/

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix (MultiCore):");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            unsigned long min_tile_size(Configuration::instance()->get_value("mc::product[DM,DM]::min_tile_size", 512));

            if (a.rows() < min_tile_size || b.columns() < min_tile_size)
            {
                return Product<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));

                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long parts_row(a.rows() / min_tile_size);
                unsigned long parts_col(b.columns() / min_tile_size);

                unsigned long i(0);
                for ( ; i + 1 < parts_row ; ++i)
                {
                    const DenseMatrixTile<DT1_> dmt_a(a, min_tile_size, a.columns(), i * min_tile_size, 0);
                    unsigned long j(0);
                    for ( ; j + 1 < parts_col ; ++j)
                    {
                        DenseMatrixTile<DT1_> dmt_r(result, min_tile_size, min_tile_size, i * min_tile_size, j * min_tile_size);
                        const DenseMatrixTile<DT1_> dmt_b(b, b.rows(), min_tile_size, 0, j * min_tile_size);

                        ThreeArgWrapper< Product<typename Tag_::DelegateTo>, DenseMatrixTile<DT1_>,
                                            const DenseMatrixTile<DT1_>, const DenseMatrixTile<DT2_> >
                                        wrapper(dmt_r, dmt_a, dmt_b);
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                    DenseMatrixTile<DT1_> dmt_r(result, min_tile_size, b.columns() % min_tile_size + min_tile_size, i * min_tile_size, j * min_tile_size);
                    const DenseMatrixTile<DT2_> dmt_b(b, b.rows(), b.columns() % min_tile_size + min_tile_size, 0, j * min_tile_size);

                    ThreeArgWrapper< Product<typename Tag_::DelegateTo>, DenseMatrixTile<DT1_>,
                                            const DenseMatrixTile<DT1_>, const DenseMatrixTile<DT2_> >
                                        wrapper(dmt_r, dmt_a, dmt_b);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
                unsigned long j(0);
                for ( ; j + 1 < parts_col ; ++j)
                {
                    DenseMatrixTile<DT1_> dmt_r(result, a.rows() % min_tile_size + min_tile_size, min_tile_size, i * min_tile_size, j * min_tile_size);
                    const DenseMatrixTile<DT1_> dmt_a(a, a.rows() % min_tile_size + min_tile_size, a.columns(), i * min_tile_size, 0);
                    const DenseMatrixTile<DT2_> dmt_b(b, b.rows(), min_tile_size, 0, j * min_tile_size);

                    ThreeArgWrapper< Product<typename Tag_::DelegateTo>, DenseMatrixTile<DT1_>,
                                            const DenseMatrixTile<DT1_>, const DenseMatrixTile<DT2_> >
                                        wrapper(dmt_r, dmt_a, dmt_b);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
                DenseMatrixTile<DT1_> dmt_r(result, a.rows() % min_tile_size + min_tile_size, b.columns() % min_tile_size + min_tile_size, i * min_tile_size, j * min_tile_size);
                const DenseMatrixTile<DT1_> dmt_a(a, a.rows() % min_tile_size + min_tile_size, a.columns(), i * min_tile_size, 0);
                const DenseMatrixTile<DT2_> dmt_b(b, b.rows(), b.columns() % min_tile_size + min_tile_size, 0, j * min_tile_size);

                ThreeArgWrapper< Product<typename Tag_::DelegateTo>, DenseMatrixTile<DT1_>,
                                        const DenseMatrixTile<DT1_>, const DenseMatrixTile<DT2_> >
                                    wrapper(dmt_r, dmt_a, dmt_b);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }

                return result;
            }

        }


        // Help function for BandedMatrix * DenseVector MultiCore
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & res_range, const BandedMatrix<DT1_> & a,
                const DenseVectorRange<DT2_> & b, unsigned long offset, unsigned long size)
        {
            int middle_index(a.rows() -1);

                for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                        vi != vi_end ; ++vi)
                {
                    // If we are below the diagonal band
                    if ((middle_index > vi.index()) && (offset + size > middle_index - vi.index()) )
                    {
                        unsigned long start(middle_index - vi.index());

                        if (offset >= start)
                        {
                            DenseVectorRange<DT1_> range_1(*vi, size, offset);
                            DenseVectorRange<DT2_> range_2(b.range(size, offset - start));
                            ScaledSum<typename Tag_::DelegateTo>::value(res_range, range_1, range_2);
                        }
                        else
                        {
                            DenseVectorRange<DT1_> range_1(*vi, offset + size - start, start);
                            DenseVectorRange<DT2_> range_2(b.range(offset + size - start, 0));
                            DenseVectorRange<DT1_> res_r(res_range, offset + size - start, start - offset);
                            ScaledSum<typename Tag_::DelegateTo>::value(res_r, range_1, range_2);
                        }
                    }
                    // If we are above or on the diagonal band
                    else if ((vi.index() >= middle_index) && (offset < vi->size() - vi.index() + middle_index))
                    {
                        unsigned long end(vi->size() - vi.index() + middle_index);

                        if (end >= offset + size)
                        {
                            DenseVectorRange<DT1_> range_1(*vi, size, offset);
                            DenseVectorRange<DT2_> range_2(b.range(size, offset + vi.index() - middle_index));
                            ScaledSum<typename Tag_::DelegateTo>::value(res_range, range_1, range_2);
                        }
                        else
                        {
                            DenseVectorRange<DT1_> range_1(*vi, end - offset, offset);
                            DenseVectorRange<DT2_> range_2(b.range(end - offset, offset + vi.index() - middle_index));
                            DenseVectorRange<DT1_> res_r(res_range, end - offset, 0);
                            ScaledSum<typename Tag_::DelegateTo>::value(res_r, range_1, range_2);
                        }
                    }
                }
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const BandedMatrix<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with DenseVector(MultiCore):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::product[BM,DV]::min-part-size", 1024));
            unsigned long overall_size(a.rows());

            DenseVector<DT1_> result(a.rows(), DT1_(0));

            if (overall_size < (min_part_size << 1))
            {
                DenseVectorRange<DT1_> _result(result.range(overall_size, 0));
                const DenseVectorRange<DT2_> _b(b.range(b.size(), 0));
                MCProduct<Tag_>::value(_result, a, _b, (unsigned long) 0, a.rows());
                return result;
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::product[BM,DVCB]::max-count", num_threads));
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16,  a.rows(), PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    DenseVectorRange<DT1_> res_range(result, part_size, offset);
                    FiveArgWrapper<MCProduct<Tag_>, DenseVectorRange<DT1_>, const BandedMatrix<DT1_>, 
                    const DenseVectorRange<DT2_>, const unsigned long, const unsigned long> mywrapper(res_range, a, b.range(b.size(), 0), offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::product[SM,DM]::min-part-size", 256));
            unsigned long overall_size(a.rows() * a.columns());

            if (overall_size < (min_part_size << 1))
            {
                return Product<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));

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
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
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
        }
    };
}
#endif
