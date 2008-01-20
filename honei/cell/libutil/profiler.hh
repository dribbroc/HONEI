/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Util C++ library. LibUtil is free software;
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

#ifndef CELL_GUARD_UTIL_PROFILER_HH
#define CELL_GUARD_UTIL_PROFILER_HH 1

#include <honei/cell/interface.hh>

#include <spu_mfcio.h>

namespace honei
{
    namespace cell
    {
        namespace intern
        {
            /**
             * \name Profiler stack
             * \{
             */

            /// Our top of the stack position
            extern unsigned profiler_top;

            /// Our address stack.
            extern LocalStoreAddress profiler_addresses[16];

            /// Our decrementer stack;
            extern unsigned profiler_decrementers[16];

            /// \}
        }

        inline void profiler_start() __attribute__((always_inline));

        inline void profiler_stop() __attribute__((always_inline));

        inline void profiler_start()
        {
#ifdef HONEI_PROFILER
            asm("brsl %0, 1f\n"
                "1: nop\n"
                : "=r" (intern::profiler_addresses[intern::profiler_top])
                : /* none */
                : /* none */
            );

            intern::profiler_decrementers[intern::profiler_top] = spu_read_decrementer();
            spu_write_decrementer(0xFFFFFFFFU);

            ++intern::profiler_top;
#endif
        }

        inline void profiler_stop()
        {
#ifdef HONEI_PROFILER
            --intern::profiler_top;

            unsigned decrementer(0xFFFFFFFFU - spu_read_decrementer());
            spu_write_decrementer(intern::profiler_decrementers[intern::profiler_top]);

            spu_write_out_intr_mbox(km_profiler);
            spu_write_out_mbox(reinterpret_cast<unsigned>(intern::profiler_addresses[intern::profiler_top]));
            spu_write_out_mbox(decrementer);
#endif
        }
    }
}

#endif
