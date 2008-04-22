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

#ifndef LIBUTIL_GUARD_TAGS_HH
#define LIBUTIL_GUARD_TAGS_HH 1

#include <honei/util/instantiation_policy.hh>

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
            tv_cell,
            tv_gpu,
            tv_fake /* used by unit tests */
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
            const static std::string name;

            /**
             * Tag-type for SSE1/2-optimised operations.
             *
             * \ingroup grptagscpusse
             */
            struct SSE :
                public InstantiationPolicy<CPU::SSE, NonCopyable>
            {
                const static TagValue tag_value = tv_cpu;
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
                const static std::string name;

                typedef tags::CPU DelegateTo;

                /**
                 * Tag-type for SSE1/2-optimised multithreaded operations.
                 *
                 * \ingroup grptagscpumulticore
                 */
                struct SSE :
                    public InstantiationPolicy<MultiCore::SSE, NonCopyable>
                {
                    const static TagValue tag_value = tv_cpu_multi_core;
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
            const static std::string name;
        };
    }

    /**
     * Output operator for tags::TagValue.
     */
    std::ostream & operator<< (std::ostream & left, tags::TagValue value);
}

#endif
