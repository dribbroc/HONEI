/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 David Gies <david-gies@gmx.de>
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

#ifndef LIBLA_GUARD_ELEMENT_PRODUCT_MC_HH
#define LIBLA_GUARD_ELEMENT_PRODUCT_MC_HH 1

#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/tags.hh>

#include <honei/la/dense_vector_range.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>

namespace honei
{
    template <typename Tag_> struct ElementProduct;

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
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of  DenseVectorContinuousBase elements (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_product[DVCB,DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::element_product[DVCB,DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT1_> range_1(a.range(part_size, offset));
                    DenseVectorRange<DT2_> range_2(b.range(part_size, offset));
                    TwoArgWrapper<ElementProduct<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVectorContinuousBase elements (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_product[SV,DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.used_elements());
            if (overall_size < 2 * min_part_size)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::element_product[SV,DVCB]::max-count", num_threads ));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list<std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                typename Vector<DT1_>::ConstElementIterator r(a.begin_non_zero_elements());

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    part_size = p->size;
                    offset = r.index();
                    r += (p->size - 1);
                    DenseVectorRange<DT1_> range(b.range(r.index() - offset + 1, offset));
                    ThreeArgWrapper<ElementProduct<Tag_>, SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(a, range, p->start);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                    ++r;
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
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

            //ThreadPool * p(ThreadPool::instance());
            PoolTask * pt[2 * a.size()-1];
            int taskcount(0);
            typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_bands());
            for (typename BandedMatrix<DT1_>::BandIterator l(a.begin_bands()),
                    l_end(a.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists())
                {
                    ++r;
                    continue;
                }

                DenseVector<DT1_> band(*l);
                TwoArgWrapper< ElementProduct<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(band, *r);
                pt[taskcount] = ThreadPool::instance()->dispatch(mywrapper);
                ++r;
                ++taskcount;
            }
            for (int j = 0 ; j < taskcount ; ++j )
            {
                pt[j]->wait_on();
                delete pt[j];
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_product[DM,DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::element_product[DM,DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                for (int i(0) ; i < a.rows() ; ++i)
                {
                    ElementProduct<Tag_>::value(a[i], b[i]);
                }
            }
            else if (overall_size < 2)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    FourArgWrapper< ElementProduct<Tag_>, DenseMatrix<DT1_>, const DenseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_product[SM,DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::element_product[SM,DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    FourArgWrapper< ElementProduct<Tag_>, SparseMatrix<DT1_>, const DenseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static void value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_product[SM,SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::element_product[SM,SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    FourArgWrapper< ElementProduct<Tag_>, SparseMatrix<DT1_>, const SparseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }
        
        template <typename DT1_, typename DT2_>
        static void value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                ElementProduct<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }
    };
}
#endif
