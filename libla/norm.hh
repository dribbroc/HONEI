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

#ifndef LIBLA_GUARD_NORM_HH
#define LIBLA_GUARD_NORM_HH 1

#include <libla/vector.hh>
#include <libla/dot_product.hh>
#include <libutil/tags.hh>

#include <cmath>

namespace honei
{
    /**
     * A VectorNormType is a template tag parameter for the VectorNorm class
     * template. It governs the type of the (mathematical) norm that shall be
     * computed.
     *
     * \ingroup grplaoperations
     */
    enum VectorNormType {
        vnt_max = 0, /// < Maximum or infinity norm.
        vnt_l_one, /// < One norm.
        vnt_l_two, /// < Two norm.
        // Extend here if higher norms are needed
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

        template <typename DT_>
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_non_zero_elements()), l_end(x.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                if (fabs(*l) > result)
                {
                    result = fabs(*l);
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_non_zero_elements()), l_end(x.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                result += fabs(*l);
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
        static DT_ value(const SparseVector<DT_> & x)
        {
            CONTEXT("When calculating norm of a SparseVector:");
            return sqrt(Norm<vnt_l_two, false>::value(x));
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_elements()), l_end(x.end_elements()) ;
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

            for (typename Vector<DT_>::ConstElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                result += exp(k * log(fabs(*l)));
            }

            return exp(log(result) / k);
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
    template <> struct Norm<vnt_l_two, false, tags::CPU::SSE>
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

        static float value(const DenseVector<float> & x)
        {
            return DotProduct<tags::CPU::SSE>::value(x, x);
        };

        static double value(const DenseVector<double> & x)
        {
            return DotProduct<tags::CPU::SSE>::value(x, x);
        };

        /// \}
    };
}

#endif
