/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_SPE_MANAGER_HH
#define LIBUTIL_GUARD_SPE_MANAGER_HH 1

#include <honei/cell/cell.hh>
#include <honei/libutil/assertion.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/instantiation_policy.hh>
#include <honei/libutil/memory_backend.hh>
#include <honei/libutil/private_implementation_pattern.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/spe.hh>

#include <string>
#include <vector>
#include <tr1/memory>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

namespace honei
{
    class SPEInstruction;
    class SPEInstructionQueue;

    /**
     * SPEManager handles all available SPEs and dispatches tasks to them.
     *
     * \ingroup grpcell
     */
    class SPEManager :
        public InstantiationPolicy<SPEManager, Singleton>,
        public PrivateImplementationPattern<SPEManager, Single>
    {
        private:
            /// \name Basic Operations
            /// \{

            /// Constructor.
            SPEManager();

            /// Destructor.
            ~SPEManager();

            /// \}

        public:
            friend class InstantiationPolicy<SPEManager, Singleton>;
            friend class SPEInstruction;
            friend class SPEInstructionQueue;

            /// \name Iterate over our SPEs
            /// \{

            typedef libwrapiter::ForwardIterator<SPEManager, SPE> Iterator;

            Iterator begin() const;
            Iterator end() const;
            unsigned int spe_count() const;

            /// \}

            /// Dispatch an SPETask to a SPE.
            void dispatch(const SPEInstruction & instruction);

            // Dispatch all SPEInstructions in the queue to one single SPE.
            void dispatch(SPEInstructionQueue & instruction_queue);

            // Wait for a SPEinstructin to be finished.
            void wait(SPEInstruction & instruction);

            // Wait for all SPEInstrictions in the queue to finish.
            void wait(SPEInstructionQueue & instruction_queue);

            // Returns the count of available spe's.
            unsigned long spe_count();
    };
}

#endif
