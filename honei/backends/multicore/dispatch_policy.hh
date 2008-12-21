/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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

#ifndef MULTICORE_GUARD_POLICY_HH
#define MULTICORE_GUARD_POLICY_HH 1

namespace honei
{
    namespace mc
    {
        // Forward declaration for friend class setting
        class ThreadPool;

        class DispatchPolicy
        {
            private:
               unsigned core_nr;
               unsigned lfd_nr;
               bool same;

               DispatchPolicy(unsigned core, unsigned lfd, unsigned s) :
                   core_nr(core),
                   lfd_nr(lfd),
                    same(s)
                {
                }

            public:

                friend class ThreadPool;

                /// Named constructors

                /// Dispatch on any core available
                static DispatchPolicy any_core()
                {
                    return DispatchPolicy(0, 0, 0);
                }

                /// Dispatch on same core as earlier task
                static DispatchPolicy same_core_as(unsigned lfd_nr)
                {
                    return DispatchPolicy(0, lfd_nr, 1);
                }

                /// Dispatch on explicit core
                static DispatchPolicy on_core(unsigned core_nr)
                {
                    return DispatchPolicy(core_nr, 0, 0);
                }

                /// Not working yet
                static DispatchPolicy not_same_core_as(unsigned lfd_nr)
                {
                    return DispatchPolicy(0, lfd_nr, 0);
                }
        };
    }
}
#endif
