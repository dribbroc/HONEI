/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_REDUCTION_HH
#define LIBLA_GUARD_REDUCTION_HH 1

#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/reduction-fwd.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/matrix_error.hh>
#include <honei/libutil/tags.hh>

#include <honei/libutil/pool_task.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libutil/wrapper.hh>
#include <honei/libutil/configuration.hh>
#include <honei/libutil/partitioner.hh>
#include <honei/libutil/benchmark_info.hh>

namespace honei
{
    template <ReductionType type_, typename Tag_ = tags::CPU> struct Reduction;

    /**
     * \brief Reduction of an entity to the sum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the sum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \sum(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum>
    {
        /**
         * \name Reductions
         * \{
         *
         * \brief Returns the scalar which is the sum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x If a is a matrix, a vector will be returned, else it will be a scalar.
         *
         */

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows());

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo VectorIterator!
            {
                result[i] = Reduction<rt_sum>::value(a[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & a)
        {
            CONTEXT("When reducing SparseMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows());

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo VectorIterator!
            {
                /// \ todo use optimized Reduction(SparseVector) instead of Reduction(Vector)
                result[i] = Reduction<rt_sum>::value(a[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & a)
        {
            CONTEXT("When reducing BandedMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows(), DT_(0));

            for (typename Matrix<DT_>::ConstElementIterator i(a.begin_elements()),
                    i_end(a.end_elements()) ; i != i_end ; ++i)
            {
                result[i.row()] += *i;
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const Vector<DT_> & vector)
        {
            CONTEXT("When reducing Vector to Scalar by sum:");

            DT_ result(0);

            for (typename Vector<DT_>::ConstElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by sum:");

            DT_ result(0);

            for (typename Vector<DT_>::ConstElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        /// \}

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a)
        {
            BenchmarkInfo result;
            DenseVector<DT1_> temp(a.columns());
            result = Reduction<rt_sum>::get_benchmark_info(temp) * a.rows();
            result.size.push_back(a.rows() * a.columns());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT1_> & a)
        {
            BenchmarkInfo result;
            for (unsigned long i(0) ; i < a.rows() ; ++i)
            {
                result = result + Reduction<rt_sum>::get_benchmark_info(a[i]);
                result.scale += a[i].used_elements();
            }
            result.size.push_back(a.rows() * a.columns());
            result.scale = (double(a.rows() * a.columns()) / result.scale);
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorBase<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = sizeof(DT1_) * a.size();
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(SparseVector<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.used_elements();
            result.load = sizeof(DT1_) * a.used_elements();
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }
    };

    /**
     * \brief Reduction of an entity to the maximum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the maximum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \max(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_max>
    {
        /**
         * \brief Returns the scalar which is the maximum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x If a is a matrix, a vector will be returned, else it will be a scalar.
         */

        /// \{
        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing DenseMatrix to Vector by max");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_max>::value(matrix[i]);
            }

            return result;
        }

                template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing SparseMatrix to Vector by max");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_max>::value(matrix[i]);
            }

            return result;
        }

                template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing BandedMatrix to Vector by max");
            /// \todo Use band interator.
            DenseVector<DT_> result(matrix.rows());

            for (typename Matrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
                    i_end(matrix.end_elements()) ; i != i_end ; ++i)
            {
                if (i.column()==0)
                    result[i.row()] = *i;
                else if (*i > result[i.row()])
                    result[i.row()] = *i;
            }
            return result;
        }

                template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & vector)
        {
            CONTEXT("When reducing DenseVectorBase to Scalar by max");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l > result)
                {
                    result = *l;
                }
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by max");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l > result)
                {
                    result = *l;
                }
            }

            return result;
        }

        /// \}
    };

    /**
     * \brief Reduction of an entity to the minimum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the minimum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \min(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_min>
    {
        /**
         * \name Reductions
         * \{
         *
         * \brief Returns the scalar which is the minimum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x Will return scalar x.
         */

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing DenseMatrix to Vector by min");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_min>::value(matrix[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing SparseMatrix to Vector by min");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_min>::value(matrix[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing BandedMatrix to Vector by max");
            /// \todo Use band interator.

            DenseVector<DT_> result(matrix.rows());

            for (typename Matrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
                    i_end(matrix.end_elements()) ; i != i_end ; ++i)
            {
                if (i.column() == 0)
                {
                    result[i.row()] = *i;
                }
                else if (*i < result[i.row()])
                {
                    result[i.row()] = *i;
                }
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & vector)
        {
            CONTEXT("When reducing DenseVector to Scalar by min");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l < result)
                {
                    result = *l;
                }
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by min");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l < result)
                {
                    result = *l;
                }
            }

            return result;
        }

        /// \}
    };

    /**
     * \brief reduction of a vector to scalar (sum)
     *
     * Reduction is the class template for the operation
     * \f[
     *     \texttt{Reduction}(x): \quad r \leftarrow \sum x_i,
     * \f]
     * which yields the reduction of the given vector x .
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum,tags::CPU::SSE>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the sum-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);

        static double value(const DenseVectorContinuousBase<double> & a);

        static DenseVector<float> value(const DenseMatrix<float> & a);

        static DenseVector<double> value(const DenseMatrix<double> & a);

        static float value(const SparseVector<float> & a);

        static double value(const SparseVector<double> & a);

        static DenseVector<float> value(const SparseMatrix<float> & a);

        static DenseVector<double> value(const SparseMatrix<double> & a);

    };
    /**
     * \brief reduction of a vector to scalar (sum)
     *
     * Reduction is the class template for the operation
     * \f[
     *     \texttt{Reduction}(x): \quad r \leftarrow \sum x_i,
     * \f]
     * which yields the reduction of the given vector x .
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the sum-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        static double value(const DenseVectorContinuousBase<double> & a);
        static DenseVector<double> value(const DenseMatrix<double> & a);
        static double value(const SparseVector<double> & a);
        static DenseVector<double> value(const SparseMatrix<double> & a);


        /// \}
    };

        /**
        * \brief reduction of a vector to scalar (min)
        *
        * Reduction is the class template for the operation
        * \f[
        *     \texttt{Reduction}(x): \quad r \leftarrow \min x_i,
        * \f]
        * which yields the reduction of the given vector x .
        *
        * \ingroup grplaoperations
        * \ingroup grplavectoroperations
        */
    template <> struct Reduction<rt_min,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the min-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        /// \}
    };

        /**
        * \brief reduction of a vector to scalar (max)
        *
        * Reduction is the class template for the operation
        * \f[
        *     \texttt{Reduction}(x): \quad r \leftarrow \max x_i,
        * \f]
        * which yields the reduction of the given vector x .
        *
        * \ingroup grplaoperations
        * \ingroup grplavectoroperations
        */
    template <> struct Reduction<rt_max,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the max-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        /// \}
    };

    template <ReductionType type_, typename Tag_> struct MCReduction;
    template <typename Tag_> struct MCReduction<rt_sum, Tag_>
    {
        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & a)
        {
            CONTEXT("When reducing BandedMatrix to DenseVector by sum (MultiCore):");
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[BM]::min-part-size", 64));
            unsigned long overall_size(a.size());
            DenseVector<DT_> result(a.size(), DT_(0));
            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                
                if (a.band_at(a.size()-1).exists())
                    result = a.band(0);

                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[BM]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    ThreeArgWrapper< Reduction<rt_sum, Tag_>, DenseVectorRange<DT_>, const BandedMatrix<DT_>, const unsigned long> mywrapper(range, a, offset);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(DenseVectorRange<DT_> & range, const BandedMatrix<DT_> & a, const unsigned long start)
        {
            unsigned long range_size(range.size());
            unsigned long startindex(a.size() - (start + range_size)), endindex((2 * a.size() - 1) - start);
            typename DenseVectorRange<DT_>::ElementIterator range_end(range.end_elements());
            typename BandedMatrix<DT_>::ConstVectorIterator vi(a.band_at(startindex)), end(a.band_at(a.size()-1));
            // Calculation for lower part
            for ( ; vi != end ; ++vi)
            {
                if (!vi.exists())
                    continue;
                DenseVectorRange<DT_> band(vi->range(range_size, start));
                typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());
                for (int i(0) ; (i < (long(startindex - vi.index()) + long(range_size-1))) ; ++i)
                {
                    ++r;
                    ++b;
                }
                for ( ; r != range_end ; ++r, ++b)
                {
                    *r += *b;
                }
            }
            // Calculation for upper part
            if (endindex > vi.index() + 1)
            {
                ++ vi;
                end = a.band_at(endindex);

                for ( ; vi != end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    DenseVectorRange<DT_> band(vi->range(range_size, start));
                    typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                    typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());        
                    unsigned long i = range_size;                    
                    if (i > (endindex - vi.index())) 
                        i = endindex - vi.index();
                    for ( ; i > 0 ;--i, ++r, ++b)
                    {
                        *r += *b;
                    }
                }
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by sum (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_sum, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_sum, Tag_>, const DenseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const DenseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & a)
        {
            CONTEXT("When reducing SparseMatrix to DenseVector by sum (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_sum, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_sum, Tag_>, const SparseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const SparseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & a)
        {
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++pri)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(a.range(part_size, offset));
                    ResultOneArgWrapper< Reduction<rt_sum, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > mywrapper(*pri, range);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & a)
        {
            CONTEXT("When reducing DenseVectorSlice to Scalar by sum (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[DVS]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[DVS]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++i)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_elements()), stop(a.begin_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_sum, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & a)
        {
            CONTEXT("When reducing SparseVector to Scalar by sum (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_sum[SV]::min-part-size", 64));
            unsigned long overall_size(a.used_elements());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_sum[SV]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++i)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_non_zero_elements()), stop(a.begin_non_zero_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_sum, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_sum, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DTP_, typename IT1_, typename IT2_>
        static void value(DTP_ result, IT1_ & x, IT2_ & x_end)
        {
            CONTEXT("When reducing iterator-based by sum:");

            for ( ; x != x_end ; ++x)
            {
                *result += *x;
            }
        }
    };

    template <typename Tag_> struct MCReduction<rt_max, Tag_>
    {
        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & a)
        {
            CONTEXT("When reducing BandedMatrix to DenseVector by max (MultiCore):");
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[BM]::min-part-size", 64));
            unsigned long overall_size(a.size());
            DenseVector<DT_> result(a.size(), DT_(0));
            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                
                if (a.band_at(a.size()-1).exists())
                    result = a.band(0);

                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[BM]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    ThreeArgWrapper< Reduction<rt_max, Tag_>, DenseVectorRange<DT_>, const BandedMatrix<DT_>, const unsigned long> mywrapper(range, a, offset);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(DenseVectorRange<DT_> & range, const BandedMatrix<DT_> & a, const unsigned long start)
        {
            unsigned long range_size(range.size());
            unsigned long startindex(a.size() - (start + range_size)), endindex((2 * a.size() - 1) - start);
            typename DenseVectorRange<DT_>::ElementIterator range_end(range.end_elements());
            typename BandedMatrix<DT_>::ConstVectorIterator vi(a.band_at(startindex)), end(a.band_at(a.size()-1));
            // Calculation for lower part
            for ( ; vi != end ; ++vi)
            {
                if (!vi.exists())
                    continue;
                DenseVectorRange<DT_> band(vi->range(range_size, start));
                typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());
                for (int i(0) ; (i < (long(startindex - vi.index()) + long(range_size-1))) ; ++i)
                {
                    ++r;
                    ++b;
                }
                for ( ; r != range_end ; ++r, ++b)
                {
                    if (*r < *b)
                        *r = *b;
                }
            }
            // Calculation for upper part
            if (endindex > vi.index() + 1)
            {
                ++ vi;
                end = a.band_at(endindex);

                for ( ; vi != end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    DenseVectorRange<DT_> band(vi->range(range_size, start));
                    typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                    typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());        
                    unsigned long i = range_size;                    
                    if (i > (endindex - vi.index())) 
                        i = endindex - vi.index();
                    for ( ; i > 0 ;--i, ++r, ++b)
                    {
                        if (*r < *b)
                            *r = *b;
                    }
                }
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by max (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_max, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_max, Tag_>, const DenseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const DenseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_max, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & a)
        {
            CONTEXT("When reducing SparseMatrix to DenseVector by max (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_max, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_max, Tag_>, const SparseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const SparseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_max, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & a)
        {
            CONTEXT("When reducing DenseVectorContinuousBase to Scalar by max (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++pri)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(a.range(part_size, offset));
                    ResultOneArgWrapper< Reduction<rt_max, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > mywrapper(*pri, range);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_max, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & a)
        {
            CONTEXT("When reducing DenseVectorSlice to Scalar by max (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[DVS]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[DVS]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++i)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_elements()), stop(a.begin_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_max, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_max, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }


        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & a)
        {
            CONTEXT("When reducing SparseVector to Scalar by max (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_max[SV]::min-part-size", 64));
            unsigned long overall_size(a.used_elements());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_max, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_max[SV]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++i)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_non_zero_elements()), stop(a.begin_non_zero_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_max, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_max, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DTP_, typename IT1_, typename IT2_>
        static void value(DTP_ result, IT1_ & x, IT2_ & x_end)
        {
            CONTEXT("When reducing iterator-based by max:");
            *result = *x;
            ++x;
            for ( ; x != x_end ; ++x)
            {
                if (*result < *x) *result = *x;
            }
        }
    };

    template <typename Tag_> struct MCReduction<rt_min, Tag_>
    {
        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & a)
        {
            CONTEXT("When reducing BandedMatrix to DenseVector by min (MultiCore):");
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[BM]::min-part-size", 64));
            unsigned long overall_size(a.size());
            DenseVector<DT_> result(a.size(), DT_(0));
            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                
                if (a.band_at(a.size()-1).exists())
                    result = a.band(0);

                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[BM]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    ThreeArgWrapper< Reduction<rt_min, Tag_>, DenseVectorRange<DT_>, const BandedMatrix<DT_>, const unsigned long> mywrapper(range, a, offset);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(DenseVectorRange<DT_> & range, const BandedMatrix<DT_> & a, const unsigned long start)
        {
            unsigned long range_size(range.size());
            unsigned long startindex(a.size() - (start + range_size)), endindex((2 * a.size() - 1) - start);
            typename DenseVectorRange<DT_>::ElementIterator range_end(range.end_elements());
            typename BandedMatrix<DT_>::ConstVectorIterator vi(a.band_at(startindex)), end(a.band_at(a.size()-1));
            // Calculation for lower part
            for ( ; vi != end ; ++vi)
            {
                if (!vi.exists())
                    continue;
                DenseVectorRange<DT_> band(vi->range(range_size, start));
                typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());
                for (int i(0) ; (i < (long(startindex - vi.index()) + long(range_size-1))) ; ++i)
                {
                    ++r;
                    ++b;
                }
                for ( ; r != range_end ; ++r, ++b)
                {
                    if (*r > *b)
                        *r = *b;
                }
            }
            // Calculation for upper part
            if (endindex > vi.index() + 1)
            {
                ++ vi;
                end = a.band_at(endindex);

                for ( ; vi != end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    DenseVectorRange<DT_> band(vi->range(range_size, start));
                    typename DenseVectorRange<DT_>::ElementIterator r(range.begin_elements());
                    typename DenseVectorRange<DT_>::ElementIterator b(band.begin_elements());        
                    unsigned long i = range_size;                    
                    if (i > (endindex - vi.index())) 
                        i = endindex - vi.index();
                    for ( ; i > 0 ;--i, ++r, ++b)
                    {
                        if (*r > *b)
                            *r = *b;
                    }
                }
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by min (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_min, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_min, Tag_>, const DenseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const DenseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_min, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & a)
        {
            CONTEXT("When reducing SparseMatrix to DenseVector by min (MultiCore):");
            DenseVector<DT_> result(a.rows());

            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                typename Vector<DT_>::ElementIterator l(result.begin_elements());
                for (unsigned long i(0) ; i < a.rows() ; ++i, ++l)
                {
                    *l = Reduction<rt_min, Tag_>::value(a[i]);
                }
            }
            else if (overall_size < 2)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size; 
                    DenseVectorRange<DT_> range(result.range(part_size, offset));
                    FourArgWrapper< Reduction<rt_min, Tag_>, const SparseMatrix<DT_>, DenseVectorRange<DT_>, unsigned long, unsigned long> mywrapper(a, range, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return result;
        }

        template <typename DT_>
        static void value(const SparseMatrix<DT_> & a, DenseVectorRange<DT_> & b, unsigned long start, unsigned long part_size)
        {
            typename Vector<DT_>::ElementIterator l(b.begin_elements());
            for (unsigned long i(start) ; i < (start + part_size) ; ++i, ++l)
            {
                *l = Reduction<rt_min, typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & a)
        {
            CONTEXT("When reducing DenseVectorContinuousBase to Scalar by min (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++pri)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT_> range(a.range(part_size, offset));
                    ResultOneArgWrapper< Reduction<rt_min, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > mywrapper(*pri, range);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_min, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & a)
        {
            CONTEXT("When reducing DenseVectorSlice to Scalar by min (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[DVS]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[DVS]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_elements()), stop(a.begin_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_min, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_min, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & a)
        {
            CONTEXT("When reducing SparseVector to Scalar by min (MultiCore):");
            DT_ result(0);
            unsigned long min_part_size(Configuration::instance()->get_value("mc::rt_min[SV]::min-part-size", 64));
            unsigned long overall_size(a.used_elements());

            if (overall_size < 2 * min_part_size)
            {
                result = Reduction<rt_min, typename Tag_::DelegateTo>::value(a);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::rt_min[SV]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                DenseVector<DT_> preresult(partitions.size(), DT_(0));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                int i(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p, ++i)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                    pri += i;
                    typename Vector<DT_>::ConstElementIterator start(a.begin_non_zero_elements()), stop(a.begin_non_zero_elements());
                    start += offset;
                    stop += (offset + part_size);
                    ThreeArgWrapper< Reduction<rt_min, Tag_>, typename Vector<DT_>::ElementIterator, typename Vector<DT_>::ConstElementIterator, typename Vector<DT_>::ConstElementIterator> mywrapper(pri, start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::get_instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            result = Reduction<rt_min, typename Tag_::DelegateTo>::value(preresult);
            }
            return result;
        }


        template <typename DTP_, typename IT1_, typename IT2_>
        static void value(DTP_ result, IT1_ & x, IT2_ & x_end)
        {
            CONTEXT("When reducing iterator-based by min:");
            *result = *x;
            ++x;
            for ( ; x != x_end ; ++x)
            {
                if (*result > *x) *result = *x;
            }
        }
    };
    template <ReductionType type_> struct Reduction<type_, tags::CPU::MultiCore> : MCReduction<type_, tags::CPU::MultiCore> {};
    template <ReductionType type_> struct Reduction<type_, tags::CPU::MultiCore::SSE> : MCReduction<type_, tags::CPU::MultiCore::SSE> {};
}
#endif
