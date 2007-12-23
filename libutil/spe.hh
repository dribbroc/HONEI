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

#ifndef LIBUTIL_GUARD_SPE_HH
#define LIBUTIL_GUARD_SPE_HH 1

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/memory_backend.hh>
#include <libutil/stringify.hh>

#include <string>
#include <tr1/memory>
#include <tr1/functional>

#include <libspe2.h>

namespace honei
{
    class SPE;
    class SPEKernel;

    /**
     * An instance of SPE encapsulates one of the system's Synergistic
     * Processing Elements.
     *
     * \ingroup grpcell
     */
    class SPE
    {
        private:
            /// Our implementation class.
            class Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

            /// \name Constructors.
            /// \{

            /// Constructor.
            SPE();

            /// \}

        public:
            friend class SPEKernel;
            friend class SPEManager;
            SPE(const SPE & other);

            /// Destructor.
            ~SPE();

            /// Return out libspe2 context.
            spe_context_ptr_t context() const;

            /// Run an SPE kernel.
            void run(const SPEKernel & kernel);

            /// Return our device id.
            DeviceId id() const;

            /// Return true if we idle.
            bool idle() const;

            /// Return our current kernel.
            SPEKernel * kernel() const;

            /// \name Mailbox functions
            /// \{

            void send_mail(unsigned int mail) const;

            /// \}

            /// \name Synchonization functions
            /// \{

            void signal() const;

            /// \}
    };
}

#endif
