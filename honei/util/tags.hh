/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtility is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBUTIL_GUARD_TAGS_HH
#define LIBUTIL_GUARD_TAGS_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/attributes.hh>

#include <ostream>

namespace honei
{
    namespace tags
    {
        /**
         * Tag values for all supported destinations.
         */
        enum TagValue
        {
            tv_cpu = 1,
            tv_cpu_multi_core,
            tv_cpu_itanium,
            tv_cell,
            tv_gpu,
            tv_gpu_cuda,
            tv_gpu_multi_core,
            tv_opencl,
            tv_opencl_cpu,
            tv_opencl_gpu,
            tv_opencl_accelerator,
            tv_fake, /* used by unit tests */
            tv_none
        };


        /**
         * Tag-type for generic/C++-based operations.
         *
         * \ingroup grptagscpu
         */
        struct CPU :
            public InstantiationPolicy<CPU, NonCopyable>
        {
            const static TagValue tag_value = tv_cpu;
            const static TagValue memory_value = tv_cpu;
            const static std::string name;

            /**
             * Tag-type for SSE-optimised operations.
             *
             * \ingroup grptagscpusse
             */
            struct SSE :
                public InstantiationPolicy<CPU::SSE, NonCopyable>
            {
                const static TagValue tag_value = tv_cpu;
                const static TagValue memory_value = tv_cpu;
                const static std::string name;
            };

            /**
             * Tag-type for Itanium-optimised operations.
             *
             * \ingroup grptagscpuitanium
             */
            struct Itanium :
                public InstantiationPolicy<CPU::Itanium, NonCopyable>
            {
                const static TagValue tag_value = tv_cpu_itanium;
                const static TagValue memory_value = tv_cpu;
                const static std::string name;
            };

            /**
             * Tag-type for multithreaded/C++-based operations.
             *
             * \ingroup grptagscpumulticore
             */
            struct MultiCore :
                public InstantiationPolicy<MultiCore, NonCopyable>
            {
                const static TagValue tag_value = tv_cpu_multi_core;
                const static TagValue memory_value = tv_cpu;
                const static std::string name;

                typedef tags::CPU DelegateTo;

                /**
                 * Tag-type for SSE-optimised multithreaded operations.
                 *
                 * \ingroup grptagscpumulticore
                 */
                struct SSE :
                    public InstantiationPolicy<MultiCore::SSE, NonCopyable>
                {
                    const static TagValue tag_value = tv_cpu_multi_core;
                    const static TagValue memory_value = tv_cpu;
                    const static std::string name;

                    typedef tags::CPU::SSE DelegateTo;
                };
            };
        };

        /**
         * Tag-type for Cell-based operations.
         *
         * \ingroup grptagscell
         */
        struct Cell :
                public InstantiationPolicy<Cell, NonCopyable>
        {
            const static TagValue tag_value = tv_cell;
            const static TagValue memory_value = tv_cpu;

            struct SPE;

            struct PPE;

            const static std::string name;
        };

        /**
         * Tag-type for GPU-base operations.
         *
         * \ingroup gtptagsgpu
         */
        struct GPU :
                public InstantiationPolicy<GPU, NonCopyable>
        {
            const static TagValue tag_value = tv_gpu;
            const static TagValue memory_value = tv_gpu;
            const static std::string name;

            /**
             * Tag-type for CUDA-optimised operations.
             *
             * \ingroup grptagsgpu
             */
            struct CUDA :
                public InstantiationPolicy<GPU::CUDA, NonCopyable>
            {
                const static TagValue tag_value = tv_gpu_cuda;
                const static TagValue memory_value = tv_gpu_cuda;
                const static std::string name;
            };

            /**
             * Tag-type for multithreaded/C++-based multi gpu operations.
             *
             * \ingroup grptagsgpumulticore
             */
            struct MultiCore :
                public InstantiationPolicy<MultiCore, NonCopyable>
            {
                const static TagValue tag_value = tv_gpu_multi_core;
                const static TagValue memory_value = tv_gpu;
                const static std::string name;

                struct CUDA :
                    public InstantiationPolicy<MultiCore::CUDA, NonCopyable>
                {
                    const static TagValue tag_value = tv_gpu_multi_core;
                    const static TagValue memory_value = tv_gpu_cuda;
                    const static std::string name;

                    typedef tags::GPU::CUDA DelegateTo;
                };
            };
        };

        /**
         * Tag-type for OpenCL-base operations.
         *
         * \ingroup gtptagsocl
         */
        struct OpenCL :
                public InstantiationPolicy<OpenCL, NonCopyable>
        {
            const static TagValue tag_value = tv_opencl;
            const static TagValue memory_value = tv_opencl;
            const static std::string name;

            /**
             * Tag-type for CPU-optimised operations.
             *
             * \ingroup grptagsocl
             */
            struct CPU :
                public InstantiationPolicy<OpenCL::CPU, NonCopyable>
            {
                const static TagValue tag_value = tv_opencl;
                const static TagValue memory_value = tv_opencl_cpu;
                const static std::string name;
            };

            /**
             * Tag-type for GPU-optimised operations.
             *
             * \ingroup grptagsocl
             */
            struct GPU :
                public InstantiationPolicy<OpenCL::GPU, NonCopyable>
            {
                const static TagValue tag_value = tv_opencl;
                const static TagValue memory_value = tv_opencl_gpu;
                const static std::string name;
            };

            /**
             * Tag-type for Accelerator-optimised operations.
             *
             * \ingroup grptagsocl
             */
            struct Accelerator :
                public InstantiationPolicy<OpenCL::Accelerator, NonCopyable>
            {
                const static TagValue tag_value = tv_opencl;
                const static TagValue memory_value = tv_opencl_accelerator;
                const static std::string name;
            };

        };

        /**
         * Tag-type for none architecture specific operations.
         *
         * \ingroup gtptagsnone
         */
        struct NONE :
                public InstantiationPolicy<NONE, NonCopyable>
        {
            const static TagValue tag_value = tv_none;
            const static TagValue memory_value = tv_none;
            const static std::string name;
        };
    }

    /**
     * Output operator for tags::TagValue.
     */
    std::ostream & operator<< (std::ostream & left, tags::TagValue value);
}

#endif
