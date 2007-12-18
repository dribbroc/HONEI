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

#include <iostream>

#include <cmath>

///\todo: Do not use define for setting size of multicore-partitions.
// For optimization purposes
#define PARTS 8

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
std::cout << "yo  " << std::endl;
            CONTEXT("When calculating norm of a DenseVectorBase:");
            return DotProduct<tags::CPU>::value(x, x);
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
std::cout << "yo  " << std::endl;
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
         * \brief Returns the L2 norm of a vector.
         *
         * \param x The entity whose norm is to be computed.
         *
         * \retval r Will create a new object of type DT_ and return it.
         */

        static float value(const DenseVectorContinuousBase<float> & x);

        static double value(const DenseVectorContinuousBase<double> & x);

        static float value(const SparseVector<float> & x);

        static double value(const SparseVector<double> & x);
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

        static float value(const SparseVector<float> & x);

        static double value(const SparseVector<double> & x);
    };

    template <> struct Norm<vnt_max, false, tags::Cell>
    {

        static float value(const DenseVectorContinuousBase<float> & a);

    };

    template <> struct Norm<vnt_l_two, true, tags::Cell>
    {

        static float value(const DenseVectorContinuousBase<float> & a);
    };


    template <> struct Norm<vnt_l_two, false, tags::Cell>
    {

        static float value(const DenseVectorContinuousBase<float> & a);

    };

    // MCNorm
    template <VectorNormType norm_type_ = vnt_l_two, bool root_ = false, typename Tag_ = tags::CPU::MultiCore> 
    struct MCNorm
    {
    };

    template <bool root_, typename Tag_> 
    struct MCNorm<vnt_max, root_, Tag_>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
std::cout << "yo  yo  yo0dvcb" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorContinuousBase (MultiCore):");
            DT_ result(0);

            unsigned long parts(PARTS);
            unsigned long div(x.size() / parts);
            if (div <= 1) 
            {
                result = Norm<vnt_max, root_, typename Tag_::DelegateTo>::value(x);
                return result;
            }
            unsigned long modulo = x.size() % parts;
            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[parts];
            DenseVector<DT_> preresult(parts, DT_(0));
            for (unsigned long i(0) ; i < (modulo) ; ++i)
            {
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                pri += i;
                DenseVectorRange<DT_> range(x.range(div+1, i * (div + 1)));
                ResultOneArgWrapper< Norm<vnt_max, root_, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > wrapper(*pri, range);
                pt[i] = p->dispatch(wrapper);
            }
            for (unsigned long i(modulo) ; i < parts ; ++i)
            {
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                pri += i;
                DenseVectorRange<DT_> range(x.range(div, modulo + div * i));
                ResultOneArgWrapper< Norm<vnt_max, root_, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > wrapper(*pri, range);
                pt[i] = p->dispatch(wrapper);
            }
            for (unsigned long i(0) ; i < parts ; ++i)
            {
                pt[i]->wait_on();
            }
            result = Norm<vnt_max, root_, typename Tag_::DelegateTo>::value(preresult);
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & x)
        {
std::cout << "yo  yo  yo0dvs" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorSlice (MultiCore):");
            // mc->sc dummy
            return Norm<vnt_max, root_, typename Tag_::DelegateTo>::value(x);
        }        

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
std::cout << "yo  yo  yo0sv" << std::endl;
            CONTEXT("When calculating norm of a SparseVector (MultiCore):");
            // mc->sc dummy
            return Norm<vnt_max, root_, typename Tag_::DelegateTo>::value(x);
        }
    };

    template <bool root_, typename Tag_> 
    struct MCNorm<vnt_l_one, root_, Tag_>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
std::cout << "yo  yo  yo1dvcb" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorContinuousBase (MultiCore):");

            DT_ result(0);

            unsigned long parts(PARTS);
            unsigned long div(x.size() / parts);
            if (div <= 1) 
            {
                result = Norm<vnt_l_one, root_, typename Tag_::DelegateTo>::value(x);
                return result;
            }
            unsigned long modulo = x.size() % parts;
            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[parts];
            DenseVector<DT_> preresult(parts, DT_(0));
            for (unsigned long i(0) ; i < (modulo) ; ++i)
            {
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                pri += i;
                DenseVectorRange<DT_> range(x.range(div+1, i * (div + 1)));
                ResultOneArgWrapper< Norm<vnt_l_one, root_, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > wrapper(*pri, range);
                pt[i] = p->dispatch(wrapper);
            }
            for (unsigned long i(modulo) ; i < parts ; ++i)
            {
                typename Vector<DT_>::ElementIterator pri(preresult.begin_elements());
                pri += i;
                DenseVectorRange<DT_> range(x.range(div, modulo + div * i));
                ResultOneArgWrapper< Norm<vnt_l_one, root_, typename Tag_::DelegateTo>, DT_, const DenseVectorRange<DT_> > wrapper(*pri, range);
                pt[i] = p->dispatch(wrapper);
            }
            for (unsigned long i(0) ; i < parts ; ++i)
            {
                pt[i]->wait_on();
            }
            result = Norm<vnt_l_one, root_, typename Tag_::DelegateTo>::value(preresult);
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & x)
        {
std::cout << "yo  yo  yo1dvs" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorSlice (MultiCore):");
            // mc->sc dummy
            return Norm<vnt_l_one, root_, typename Tag_::DelegateTo>::value(x);
        }        

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
std::cout << "yo  yo  yo1sv" << std::endl;
            CONTEXT("When calculating norm of a SparseVector (MultiCore):");
            // mc->sc dummy
            return Norm<vnt_l_one, root_, typename Tag_::DelegateTo>::value(x);
        }
    };

    template <typename Tag_>
    struct MCNorm<vnt_l_two, false, Tag_>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
std::cout << "yo  yo  yo2dv" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorContinuousBase (MultiCore):");
            return DotProduct<Tag_>::value(x, x);
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & x)
        {
std::cout << "yo  yo  yo2dvb" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorSlice (MultiCore):");
            DenseVector<DT_> y(x.copy());
            return DotProduct<Tag_>::value(x, x);
        }
        
        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
std::cout << "yo  yo  yo2sv" << std::endl;
            CONTEXT("When calculating norm of a SparseVecto (Multicore):");
            return DotProduct<tags::CPU>::value(x, x);
        }
    };

    template <typename Tag_>
    struct MCNorm<vnt_l_two, true, Tag_>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & x)
        {
std::cout << "yo  yo  yo2tdvb" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorBase (MultiCore):");
            return sqrt(Norm<vnt_l_two, false, Tag_>::value(x));
        }

        template <typename DT_>
        static DT_ value(const DenseVectorSlice<DT_> & x)
        {
std::cout << "yo  yo  yo2tdvb" << std::endl;
            CONTEXT("When calculating norm of a DenseVectorBase (MultiCore):");
            return sqrt(Norm<vnt_l_two, false, Tag_>::value(x));
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & x)
        {
std::cout << "yo  yo  yo2tsv" << std::endl;
            CONTEXT("When calculating norm of a SparseVector (MultiCore):");
            return sqrt(Norm<vnt_l_two, false, Tag_>::value(x));
        }
    };

    template <VectorNormType norm_type_, bool root_> struct Norm<norm_type_, root_, tags::CPU::MultiCore> : MCNorm<norm_type_, root_, tags::CPU::MultiCore> {};
    template <VectorNormType norm_type_, bool root_> struct Norm<norm_type_, root_, tags::CPU::MultiCore::SSE> : MCNorm<norm_type_, root_, tags::CPU::MultiCore::SSE> {};
}



#endif
