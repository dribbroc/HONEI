/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
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

#pragma once
#ifndef LIBLA_GUARD_NORM_HH
#define LIBLA_GUARD_NORM_HH 1

#include <honei/la/dot_product.hh>
#include <honei/la/norm-fwd.hh>
#include <honei/util/tags.hh>
#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/util/configuration.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>

#include <cmath>


namespace honei
{
    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <VectorNormType norm_type_ = vnt_l_two, bool root_ = false, typename Tag_ = tags::CPU> struct Norm
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the given norm of a vector.
         *
         * \param x The vector whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        /*template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating norm of a DenseVectorBase:");
            DT_ result(0);

            result = Norm<norm_type_, root_>::value(x);

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseVector:");
            DT_ result(0);

            result = Norm<norm_type_, root_>::value(x);

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseMatrix<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseMatrix:");
            DT_ result(0);

            result = Norm<norm_type_, root_>::value(x);

            return result;
        }*/

        /// \}
    };

    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_{\infty},
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <bool root_> struct Norm<vnt_max, root_, tags::CPU>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L-Infinity (maximum) norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating norm of a DenseVectorBase:");

            DT_ result(0);

            for (typename DenseVectorBase<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
                    l != l_end ; ++l)
            {
                if (fabs(*l) > result)
                {
                    result = fabs(*l);
                }
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseVector:");

            DT_ result(0);

            for (typename SparseVector<DT_>::NonZeroConstElementIterator l(x.begin_non_zero_elements()), l_end(x.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                if (fabs(*l) > result)
                {
                    result = fabs(*l);
                }
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseMatrix<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseMatrix:");

            DT_ result(0);

            for (unsigned long row(0) ; row < x.rows() ; ++x)
            {
                for (typename SparseVector<DT_>::NonZeroConstElementIterator l(x[row].begin_non_zero_elements()), l_end(x[row].end_non_zero_elements()) ;
                        l != l_end ; ++l)
                {
                    if (fabs(*l) > result)
                    {
                        result = fabs(*l);
                    }
                }
            }

            return result;
        }

        /// \}
    };

    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_1,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <bool root_> struct Norm<vnt_l_one, root_, tags::CPU>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L1 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        template <typename DT_>
            static DT_ value(const DenseVectorBase<DT_> & x)
            {
                CONTEXT("When calculating norm of a DenseVectorBase:");

                DT_ result(0);

                for (typename DenseVectorBase<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
                        l != l_end ; ++l)
                {
                    result += fabs(*l);
                }

                return result;
            }

        template <typename DT_>
            static DT_ value(const SparseVector<DT_> & x)
            {
                CONTEXT("When calculating norm of a SparseVector:");
                DT_ result(0);

                for (typename SparseVector<DT_>::NonZeroConstElementIterator l(x.begin_non_zero_elements()), l_end(x.end_non_zero_elements()) ;
                        l != l_end ; ++l)
                {
                    result += fabs(*l);
                }

                return result;
            }

        template <typename DT_>
            static DT_ value(const SparseMatrix<DT_> & x)
            {
                CONTEXT("When calculating norm of a SparseMatrix:");
                DT_ result(0);

                for (unsigned long row(0) ; row < x.rows() ; ++row)
                {
                    for (typename SparseVector<DT_>::NonZeroConstElementIterator l(x[row].begin_non_zero_elements()), l_end(x[row].end_non_zero_elements()) ;
                            l != l_end ; ++l)
                    {
                        result += fabs(*l);
                    }
                }

                return result;
            }

        /// \}
    };

    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Norm<vnt_l_two, false, tags::CPU>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating norm of a DenseVectorBase:");
            return DotProduct<tags::CPU>::value(x, x);
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseVector:");
            return DotProduct<tags::CPU>::value(x, x);
        }

        template <typename DT_>
        static DT_ value(const SparseMatrix<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseMatrix:");
            DT_ result(0);
            for (unsigned long row(0) ; row < x.rows() ; ++row)
            {
                result += DotProduct<tags::CPU>::value(x[row], x[row]);
            }
            return result;
        }

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return MPIOps<tags::CPU>::norm_l2_false(x);
        }

        template <typename DT_>
        static BenchmarkInfo get_benchmark_info(const DenseVectorContinuousBase<DT_> & x)
        {
            BenchmarkInfo result;
            result.flops = 2 * x.size();
            result.load = x.size() * sizeof(DT_);
            result.store = sizeof(DT_);
            return result;
        }

        /// \}
    };

    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Norm<vnt_l_two, true, tags::CPU>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the square root of the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating norm of a DenseVectorBase:");
            return sqrt(Norm<vnt_l_two, false>::value(x));
        }

        template <typename DT_>
        static DT_ value(const SparseMatrix<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseMatrix:");
            return sqrt(Norm<vnt_l_two, false>::value(x));
        }

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return sqrt(MPIOps<tags::CPU>::norm_l2_false(x));
        }

        template <typename DT_>
        static BenchmarkInfo get_benchmark_info(const DenseVectorContinuousBase<DT_> & x)
        {
            BenchmarkInfo result;
            result.flops = 2 * x.size() + 1;
            result.load = x.size() * sizeof(DT_);
            result.store = sizeof(DT_);
            return result;
        }

        /// \}
    };

    /**
     * \brief Norm of an entity.
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_n,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <VectorNormType norm_type_> struct Norm<norm_type_, true, tags::CPU>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the square root of the given norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        template <typename DT_>
        static DT_ value (const DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating norm of a DenseVectorBase:");
            DT_ result(0);
            unsigned int k(static_cast<unsigned int>(norm_type_));

            for (typename DenseVectorBase<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
                    l != l_end ; ++l)
            {
                if (*l != static_cast<DT_>(0))
                {
                    result += exp(k * log(fabs(*l)));
                }
            }

            return exp(log(result) / k);
        }

        template <typename DT_>
        static DT_ value (const SparseVector<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseVector:");
            DT_ result(0);
            unsigned int k(static_cast<unsigned int>(norm_type_));

            for (typename SparseVector<DT_>::NonZeroConstElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                result += exp(k * log(fabs(*l)));
            }

            return exp(log(result) / k);
        }

        /// \}
    };

    template <> struct Norm<vnt_l_two, false, tags::CPU::Generic>
    {
        template <typename DT_>
        static inline DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
            BENCHADD((Norm<vnt_l_two, false, tags::CPU>::get_benchmark_info(x)));
            const DT_ * xe(x.elements());
            const unsigned long size(x.size());
            DT_ result(0);
            for (unsigned long i(0) ; i < size ; ++i)
            {
                result += xe[i]*xe[i];
            }

            return result;
        }
    };

    template <> struct Norm<vnt_l_two, true, tags::CPU::Generic>
    {
        template <typename DT_>
        static inline DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
            BENCHADD((Norm<vnt_l_two, true, tags::CPU>::get_benchmark_info(x)));
            const DT_ * xe(x.elements());
            const unsigned long size(x.size());
            DT_ result(0);
            for (unsigned long i(0) ; i < size ; ++i)
            {
                result += xe[i]*xe[i];
            }

            return sqrt(result);
        }
    };

    /**
     * \brief Norm of an entity (CUDA implementation).
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Norm<vnt_l_two, false, tags::GPU::CUDA>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);
        static double value(const DenseVectorContinuousBase<double> & x);

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return MPIOps<tags::GPU::CUDA>::norm_l2_false(x);
        }
    };
    template <> struct Norm<vnt_l_two, true, tags::GPU::CUDA>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the square root of the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);
        static double value(const DenseVectorContinuousBase<double> & x);

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return sqrt(MPIOps<tags::GPU::CUDA>::norm_l2_false(x));
        }
    };

    template <> struct Norm<vnt_l_two, false, tags::GPU::MultiCore::CUDA>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);
        static double value(const DenseVectorContinuousBase<double> & x);
    };
    template <> struct Norm<vnt_l_two, true, tags::GPU::MultiCore::CUDA>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the square root of the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);
        static double value(const DenseVectorContinuousBase<double> & x);
    };

    /**
     * \brief Norm of an entity (SSE Implementaion).
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Norm<vnt_l_two, false, tags::CPU::SSE>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);

        static double value(const DenseVectorContinuousBase<double> & x);

        static float value(const DenseMatrix<float> & x);

        static double value(const DenseMatrix<double> & x);

        static float value(const SparseVector<float> & x);

        static double value(const SparseVector<double> & x);

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return MPIOps<tags::CPU::SSE>::norm_l2_false(x);
        }
    };

    template <> struct Norm<vnt_l_two, true, tags::CPU::SSE>
    {
        /**
         * \name Norms
         * \{
         *
         * \brief Returns the square root of the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);

        static double value(const DenseVectorContinuousBase<double> & x);

        static float value(const DenseMatrix<float> & x);

        static double value(const DenseMatrix<double> & x);

        static float value(const SparseVector<float> & x);

        static double value(const SparseVector<double> & x);

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x)
        {
            return sqrt(MPIOps<tags::CPU::SSE>::norm_l2_false(x));
        }
    };

    template <> struct Norm<vnt_l_two, false, tags::OpenCL::CPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x);
    };

    template <> struct Norm<vnt_l_two, true, tags::OpenCL::CPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x);
    };

    template <> struct Norm<vnt_l_two, false, tags::OpenCL::GPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x);
    };

    template <> struct Norm<vnt_l_two, true, tags::OpenCL::GPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x);
    };

     /**
     * \brief Norm of an entity (Cell-implementation).
     *
     * Norm is the class template for the operation
     * \f[
     *     \texttt{Norm}(x): \quad r \leftarrow ||(x_{0}, \dots, x_{size - 1}||_2,
     * \f]
     * which yields c, the norm of entity x.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */

    template <> struct Norm<vnt_max, false, tags::Cell>
    {
        static float value(const DenseVectorContinuousBase<float> & a);
        static double value(const DenseVectorContinuousBase<double> & a)
        {
            /// \todo Remove CPU dummy.
            return Norm<vnt_max, false, tags::CPU>::value(a);
        }
    };

    template <> struct Norm<vnt_l_one, false, tags::Cell>
    {
        static float value(const DenseVectorContinuousBase<float> & a);
        static double value(const DenseVectorContinuousBase<double> & a);
    };


    template <> struct Norm<vnt_l_two, true, tags::Cell>
    {
        static float value(const DenseVectorContinuousBase<float> & a);
        static double value(const DenseVectorContinuousBase<double> & a);
    };


    template <> struct Norm<vnt_l_two, false, tags::Cell>
    {
        static float value(const DenseVectorContinuousBase<float> & a);
        static double value(const DenseVectorContinuousBase<double> & a);
    };
    // end of the template-Cell-part

    template <> struct Norm<vnt_l_two, false, tags::CPU::MultiCore>
    {
            template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::DelegateTo> >::op(result, x, min_part_size, max_count);

                return result;
            }

            template <typename DT_>
            static inline DT_ value(const DenseVectorMPI<DT_> & x)
            {
                return MPIOps<tags::CPU::MultiCore>::norm_l2_false(x);
            }
    };

    template <> struct Norm<vnt_l_two, true, tags::CPU::MultiCore>
    {
            template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::DelegateTo> >::op(result, x, min_part_size, max_count);
                result = sqrt(result);

                return result;
            }

            template <typename DT_>
            static inline DT_ value(const DenseVectorMPI<DT_> & x)
            {
                return sqrt(MPIOps<tags::CPU::MultiCore>::norm_l2_false(x));
            }
    };

    template <> struct Norm<vnt_l_two, false, tags::CPU::MultiCore::Generic>
    {
            template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::Generic::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::Generic::DelegateTo> >::op(result, x, min_part_size, max_count);

                return result;
            }
    };

    template <> struct Norm<vnt_l_two, true, tags::CPU::MultiCore::Generic>
    {
            template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::Generic::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::Generic::DelegateTo> >::op(result, x, min_part_size, max_count);
                result = sqrt(result);

                return result;
            }
    };

    template <> struct Norm<vnt_l_two, false, tags::CPU::MultiCore::SSE>
    {
        template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::SSE::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::SSE::DelegateTo> >::op(result, x, min_part_size, max_count);

                return result;
            }

            template <typename DT_>
            static inline DT_ value(const DenseVectorMPI<DT_> & x)
            {
                return MPIOps<tags::CPU::MultiCore::SSE>::norm_l2_false(x);
            }
    };

    template <> struct Norm<vnt_l_two, true, tags::CPU::MultiCore::SSE>
    {
            template <typename DT1_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + tags::CPU::MultiCore::SSE::name);


                DT1_ result(0);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                mc::Operation<honei::Norm<vnt_l_two, false, typename tags::CPU::MultiCore::SSE::DelegateTo> >::op(result, x, min_part_size, max_count);
                result = sqrt(result);

                return result;
            }

            template <typename DT_>
            static inline DT_ value(const DenseVectorMPI<DT_> & x)
            {
                return sqrt(MPIOps<tags::CPU::MultiCore::SSE>::norm_l2_false(x));
            }
    };
}
#endif
