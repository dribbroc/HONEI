/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef CELL_GUARD_UTIL_DEBUG_HH
#define CELL_GUARD_UTIL_DEBUG_HH 1

#include <spu_mfcio.h>

namespace honei
{
    namespace cell
    {
        inline void debug_get(const EffectiveAddress & ea, const LocalStoreAddress & lsa,
                const unsigned & size) __attribute__((always_inline));

        inline void debug_put(const EffectiveAddress & ea, const LocalStoreAddress & lsa,
                const unsigned & size) __attribute__((always_inline));

        inline void debug_acquire(const LocalStoreAddress & lsa) __attribute__((always_inline));

        inline void debug_release(const LocalStoreAddress & lsa) __attribute__((always_inline));

        inline void debug_dump() __attribute__((always_inline));

        inline void debug_get(const EffectiveAddress & ea, const LocalStoreAddress & lsa,
                const unsigned & size)
        {
            spu_write_out_intr_mbox(km_debug_get);
            spu_write_out_mbox(ea >> 32);
            spu_write_out_mbox(ea & 0xFFFFFFFF);
            spu_write_out_mbox(reinterpret_cast<unsigned int>(lsa));
            spu_write_out_mbox(size);
        }

        inline void debug_put(const EffectiveAddress & ea, const LocalStoreAddress & lsa,
                const unsigned & size)
        {
            spu_write_out_intr_mbox(km_debug_put);
            spu_write_out_mbox(ea >> 32);
            spu_write_out_mbox(ea & 0xFFFFFFFF);
            spu_write_out_mbox(reinterpret_cast<unsigned int>(lsa));
            spu_write_out_mbox(size);
        }

        inline void debug_acquire(const LocalStoreAddress & lsa)
        {
            spu_write_out_intr_mbox(km_debug_acquire_block);
            spu_write_out_mbox(reinterpret_cast<unsigned int>(lsa));
        }

        inline void debug_release(const LocalStoreAddress & lsa)
        {
            spu_write_out_intr_mbox(km_debug_release_block);
            spu_write_out_mbox(reinterpret_cast<unsigned int>(lsa));
        }

        inline void debug_dump()
        {
            spu_write_out_intr_mbox(km_debug_dump);
            spu_read_signal1();
        }
    }
}

#endif
