/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBUTIL_GUARD_TYPE_TRAITS_HH
#define LIBUTIL_GUARD_TYPE_TRAITS_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/attributes.hh>
#include <honei/util/memory_pool.hh>
#include <honei/util/memory_arbiter.hh>

#include <cstdlib>
#include <new>

namespace honei
{
    namespace intern
    {
        /**
         * Type traits for plain old data types.
         */
        template <typename DT_> struct PODTraits :
            public InstantiationPolicy<PODTraits<DT_>, NonInstantiable>
        {
            /**
             * Allocate memory that suffices to store count instances of
             * DT_. Memory needs to be aligned at 16 byte boundaries.
             *
             * \param count Count of instances that shall fit into the allocated
             *              memory.
             *
             * \note PODTraits::allocate always allocates a multiple of 16
             *       bytes.
             *
             * \todo Exceptions.
             */
            static inline DT_ * allocate(std::size_t count)
            {
                //void * result(0);

                //if (0 != posix_memalign(&result, 16, multiple_of_sixteen(sizeof(DT_) * count)))
                //    throw std::bad_alloc();
                void * result(MemoryPool<tags::CPU>::instance()->alloc(count * sizeof(DT_)));
                MemoryArbiter::instance()->register_address(result);

                return reinterpret_cast<DT_ *>(result);
            }

            static inline DT_ * reallocate(void * address, std::size_t count)
            {
                void * result(MemoryPool<tags::CPU>::instance()->realloc(address, count * sizeof(DT_)));
                MemoryArbiter::instance()->remove_address(address);
                MemoryArbiter::instance()->register_address(result);

                return reinterpret_cast<DT_ *>(result);
            }

            /**
             * Create count instances of DT_ at location by copying proto.
             *
             * \param location Memory location where instances of DT_ shall
             *                 be created.
             * \param count Count of instances that shall be created.
             * \param proto Prototype for the instances.
             *
             * \note Empty for plain old data types (POD).
             */
            static inline void create(DT_ * /*location*/, std::size_t /*count*/, HONEI_UNUSED const DT_ & proto = DT_(0))
            {
            }

            /**
             * Copy count instances of DT_ from source to dest.
             *
             * \param source Memory location whence to copy.
             * \param dest Memory location whither to copy.
             * \param count Count on instances that shall be copied.
             */
            static void copy(const DT_ * source, DT_ * dest, std::size_t count);

            /**
             * Fill count instances of DT_ in dest with value.
             *
             * \param dest Memory location whither to copy.
             * \param count Count of instances that shall be filled.
             * \param proto Prototype for the instances.
             */
            static void fill(DT_ * dest, std::size_t count, const DT_ & proto = DT_(0));

            /**
             * Destroy count instances of DT_ at location.
             *
             * \param location Memory location where instances of DT_
             *                 shall be destroyed.
             * \param count Count of instances that shall be destroyed.
             *
             * \note Empty for plain old data types (POD).
             */
            static inline void destroy(DT_ * /*location*/, std::size_t /*count*/)
            {
            }

            /**
             * Free memory that was previously allocated by PODTraits::allocate.
             *
             * \param location Memory location whose memory shall be freed.
             * \param count Count of instances for which memory was allocated.
             */
            static inline void free(DT_ * location, std::size_t /*count*/)
            {
                //::free(location);
                MemoryArbiter::instance()->remove_address(location);
                MemoryPool<tags::CPU>::instance()->free((void *)location);
            }
        };

        /**
         * Type traits for data types with either a non-trivial constructor,
         * assignment operator or destructor.
         */
        template <typename DT_> struct DefaultTraits :
            public InstantiationPolicy<DefaultTraits<DT_>, NonInstantiable>
        {
            /**
             * Allocate memory that suffices to store count instances of
             * DT_.
             *
             * \param count Count of instances that shall fit into the allocated
             *              memory.
             *
             * \todo Exceptions.
             */
            static inline DT_ * allocate(std::size_t count)
            {
                void * result(malloc(sizeof(DT_) * count));

                if (0 == result)
                    throw std::bad_alloc();

                return reinterpret_cast<DT_ *>(result);
            }

            /**
             * Create count instances of DT_ at location by copying proto.
             *
             * \param location Memory location where instances of DT_ shall
             *                 be created.
             * \param count Count of instances that shall be created.
             * \param proto Prototype for the instances.
             *
             * \note Empty for plain old data types (POD).
             */
            static inline void create(DT_ * location, std::size_t count, const DT_ & proto = DT_(0))
            {
                for (DT_ * i(location), * i_end(location + count) ; i != i_end ; ++i)
                {
                    new (i) DT_(proto);
                }
            }

            /**
             * Copy count instance of DT_ from source to dest.
             *
             * \param source Memory location whence to copy.
             * \param dest Memory location whither to copy.
             * \param count Count on instances that shall be copied.
             */
            static inline void copy(const DT_ * source, DT_ * dest, std::size_t count)
            {
                const DT_ * s(source);
                for (DT_ * d(dest), * d_end(dest + count) ; d != d_end ; ++d, ++s)
                {
                    *d = *s;
                }
            }

            /**
             * Fill count instances of DT_ in dest with value.
             *
             * \param dest Memory location whither to copy.
             * \param count Count of instances that shall be filled.
             * \param proto Prototype for the instances.
             */
            static void fill(DT_ * dest, std::size_t count, const DT_ & proto = DT_(0))
            {
                for (DT_ * d(dest), * d_end(dest + count) ; d != d_end ; ++d)
                {
                    *d = proto;
                }
            }

            /**
             * Destroy count instances of DT_ at location.
             *
             * \param location Memory location where instances of DT_
             *                 shall be destroyed.
             * \param count Count of instances that shall be destroyed.
             *
             * \note Empty fo plain old data types (POD).
             */
            static inline void destroy(DT_ * location, std::size_t count)
            {
                for (DT_ * i(location), * i_end(location + count) ; i != i_end ; ++i)
                {
                    i->~DT_();
                }
            }

            /**
             * Free memory that was previously allocated by DefaultTraits::allocate.
             *
             * \param location Memory location whose memory shall be freed.
             * \param count Count of instances for which memory was allocated.
             */
            static inline void free(DT_ * location, std::size_t /*count*/)
            {
                ::free(location);
            }
        };

        /**
         * Conversion for data types with either a non-trivial constructor,
         * assignment operator or destructor.
         */
        template <typename DT_> struct DefaultConversionTraits
        {
            template <typename CT_> void convert(CT_ * copy, const DT_ * orig, std::size_t count)
            {
                const DT_ * s(orig);
                for (DT_ * d(copy), * d_end(copy + count) ; d != d_end ; ++d, ++s)
                {
                    *d = *s;
                }
            }
        };

        /**
         * Conversion for POD.
         */
        template <typename DT_> struct PODConversionTraits
        {
        };
        template <> struct PODConversionTraits<float>
        {
            static void convert(double * copy, const float * orig, std::size_t count);
        };
        template <> struct PODConversionTraits<double>
        {
            static void convert(float * copy, const double * orig, std::size_t count);
        };
    }

    /**
     * Our type traits include functions to allocate, create, destroy and free
     * instances of type DT_.
     */

    /// \{

    template <typename DT_> struct TypeTraits :
        public intern::DefaultTraits<DT_>,
        public intern::DefaultConversionTraits<DT_>
    {
    };

    template <> struct TypeTraits<float> :
        public intern::PODTraits<float>,
        public intern::PODConversionTraits<float>
    {
    };

    template <> struct TypeTraits<double> :
        public intern::PODTraits<double>,
        public intern::PODConversionTraits<double>
    {
    };

    template <> struct TypeTraits<unsigned long> :
        public intern::PODTraits<unsigned long>,
        public intern::DefaultConversionTraits<unsigned long>
    {
    };

    template <> struct TypeTraits<long> :
        public intern::PODTraits<long>,
        public intern::DefaultConversionTraits<long>
    {
    };

    template <> struct TypeTraits<bool> :
        public intern::PODTraits<bool>,
        public intern::DefaultConversionTraits<bool>
    {
    };
    /// \}

}

#endif
