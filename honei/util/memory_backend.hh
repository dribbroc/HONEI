/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef UTIL_GUARD_MEMORY_BACKEND_HH
#define UTIL_GUARD_MEMORY_BACKEND_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/memory_backend_base.hh>
#include <honei/util/tags.hh>
#include <honei/util/exception.hh>
#include <honei/util/type_traits.hh>

#include <cstring>
#include <map>
#include <string>

namespace honei
{
    /**
     * DeviceId uniquely identifies one of the available remote memory devices that a backend
     * can interface. Possible are e.g. several Cell BE SPE Local Stores or several video memories.
     */
    typedef unsigned long DeviceId;

    /**
     * If there is only one device supported by a MemoryBackend, it must not accept requests for
     * other devices except default_device.
     */
    const DeviceId default_device(~0x0L);

    /**
     * MemoryBackend handles local copies of memory chunks in a specific device memory.
     *
     * \ingroup grpmemorymanager
     */
    template<typename Tag_>
    class MemoryBackend :
        public MemoryBackendBase
    {
    };

    template<>
    class MemoryBackend<tags::CPU> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::CPU>, Singleton>
    {
    public:
        virtual void * upload(void * /*memid*/, void * address, unsigned long /*bytes*/)
        {
            CONTEXT("When uploading data (CPU):");

            return address;
        }

        virtual void download(void * /*memid*/, void * /*address*/, unsigned long /*bytes*/)
        {
            CONTEXT("When downloading data (CPU):");
        }

        virtual void * alloc(void * /*memid*/, void * address, unsigned long /*bytes*/)
        {
            CONTEXT("When allocating data (CPU):");

            return address;
        }

        virtual void free(void * /*memid*/)
        {
           CONTEXT("When freeing data (CPU):");
        }

        virtual void copy(void * /*src_id*/, void * src_address, void * /*dest_id*/,
                void * dest_address, unsigned long bytes)
        {
           CONTEXT("When copying data (CPU):");
           std::memcpy((char *)dest_address, (char *)src_address, bytes);
        }

        virtual void convert_float_double(void * /*src_id*/, void * src_address, void * /*dest_id*/,
                void * dest_address, unsigned long bytes)
        {
           CONTEXT("When converting float->double data (CPU):");
           TypeTraits<float>::convert((double*) dest_address, (float*)src_address, bytes/sizeof(float));
        }

        virtual void convert_double_float(void * /*src_id*/, void * src_address, void * /*dest_id*/,
                void * dest_address, unsigned long bytes)
        {
           CONTEXT("When converting double->float data (CPU):");
           TypeTraits<double>::convert((float*) dest_address, (double*)src_address, bytes/sizeof(double));
        }

        virtual void fill(void * /*memid*/, void * address, unsigned long bytes, float proto)
        {
            CONTEXT("When filling data (CPU):");
            TypeTraits<float>::fill((float *)address, bytes / 4, proto);
        }

        virtual bool knows(void * /*memid*/, void * /*address*/)
        {
            return true;
        }

    };

    template<>
    class MemoryBackend<tags::GPU::CUDA> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>,
        public PrivateImplementationPattern<MemoryBackend<tags::GPU::CUDA>, Single>
    {
        public:
            /// Constructor.
            MemoryBackend<tags::GPU::CUDA>();

            /// Destructor.
            virtual ~MemoryBackend<tags::GPU::CUDA>();

            virtual void * upload(void * memid, void * address, unsigned long bytes);

            virtual void download(void * memid, void * address, unsigned long bytes);

            virtual void * alloc(void * memid, void * address, unsigned long bytes);

            virtual void free(void * memid);

            virtual void copy(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

            virtual bool knows(void * memid, void * address);
    };

    template<>
    class MemoryBackend<tags::OpenCL::CPU> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::OpenCL::CPU>, Singleton>,
        public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::CPU>, Single>
    {
        public:
            /// Constructor.
            MemoryBackend<tags::OpenCL::CPU>();

            /// Destructor.
            virtual ~MemoryBackend<tags::OpenCL::CPU>();

            virtual void * upload(void * memid, void * address, unsigned long bytes);

            virtual void download(void * memid, void * address, unsigned long bytes);

            virtual void * alloc(void * memid, void * address, unsigned long bytes);

            virtual void free(void * memid);

            virtual void copy(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

            virtual bool knows(void * memid, void * address);
    };

    template<>
    class MemoryBackend<tags::OpenCL::GPU> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::OpenCL::GPU>, Singleton>,
        public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::GPU>, Single>
    {
        public:
            /// Constructor.
            MemoryBackend<tags::OpenCL::GPU>();

            /// Destructor.
            virtual ~MemoryBackend<tags::OpenCL::GPU>();

            virtual void * upload(void * memid, void * address, unsigned long bytes);

            virtual void download(void * memid, void * address, unsigned long bytes);

            virtual void * alloc(void * memid, void * address, unsigned long bytes);

            virtual void free(void * memid);

            virtual void copy(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

            virtual bool knows(void * memid, void * address);
    };

    template<>
    class MemoryBackend<tags::OpenCL::Accelerator> :
        public MemoryBackendBase,
        public InstantiationPolicy<MemoryBackend<tags::OpenCL::Accelerator>, Singleton>,
        public PrivateImplementationPattern<MemoryBackend<tags::OpenCL::Accelerator>, Single>
    {
        public:
            /// Constructor.
            MemoryBackend<tags::OpenCL::Accelerator>();

            /// Destructor.
            virtual ~MemoryBackend<tags::OpenCL::Accelerator>();

            virtual void * upload(void * memid, void * address, unsigned long bytes);

            virtual void download(void * memid, void * address, unsigned long bytes);

            virtual void * alloc(void * memid, void * address, unsigned long bytes);

            virtual void free(void * memid);

            virtual void copy(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_float_double(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void convert_double_float(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes);

            virtual void fill(void * memid, void * address, unsigned long bytes, float proto);

            virtual bool knows(void * memid, void * address);
    };
}
#endif
