/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
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

#ifndef CELL_LIBUTIL_GUARD_ALLOCATOR_HH
#define CELL_LIBUTIL_GUARD_ALLOCATOR_HH 1

#include <cell/cell.hh>
#include <cell/libutil/debug.hh>

namespace honei
{
    namespace cell
    {
        class Allocation
        {
            private:
                unsigned used : 1;
                unsigned user_data : 31;
            public:
                friend void init_alloc(const Environment & env);
                friend Allocation * acquire_block();
                friend void release_block(Allocation & i);
                friend void release_all_blocks();
                friend unsigned get_user_data(const Allocation & i);
                friend void set_user_data(Allocation & i, unsigned data);

                void * address;
        };

        namespace intern
        {
            extern Allocation allocations[16]; /// \todo remove hardcoded numbers
        }

        inline void init_alloc(const Environment & env)
        {
            char * last_address(const_cast<char *>(reinterpret_cast<const char *>(env.begin)));

            for (Allocation * i(intern::allocations), * i_end(intern::allocations + 16) ;
                    i != i_end ; ++i)
            {
                i->used = 0;
                i->user_data = 0;
                i->address = last_address;
                last_address += Traits<honei::tags::Cell::SPE>::mfc_max_transfer_size;
            }
        }

        inline Allocation * acquire_block()
        {
            Allocation * result(0);

            for (Allocation * i(intern::allocations), * i_end(intern::allocations + 16) ; i != i_end ; ++i)
            {
                if (i->used == 1)
                    continue;

                result = i;
                result->used = 1;
                debug_acquire(result->address);
                break;
            }

            return result;
        }

        inline void release_all_blocks()
        {
            for (Allocation * i(intern::allocations), * i_end(intern::allocations + 16) ; i != i_end ; ++i)
            {
                if (i->used)
                    debug_release(i->address);

                i->user_data = 0;
                i->used = 0;
            }
        }

        inline void release_block(Allocation & i)
        {
            debug_release(i.address);
            i.user_data = 0;
            i.used = 0;
        }

        inline unsigned get_user_data(const Allocation & i)
        {
            return i.user_data;
        }

        inline void set_user_data(Allocation & i, unsigned data)
        {
            i.user_data = data & ~(1 << 31);
        }
    }
}

#endif
