/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <cell/cell.hh>
#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/memory_backend.hh>
#include <libutil/stringify.hh>
#include <libutil/spe.hh>

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
    class SPEInstructionStream;

    /**
     * SPEManager handles all available SPEs and dispatches tasks to them.
     *
     * \ingroup grpcell
     */
    class SPEManager
    {
        private:
            /// Our implementation class.
            class Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// \name Basic Operations
            /// \{

            /// Constructor.
            SPEManager();

            /// Destructor.
            ~SPEManager();

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            SPEManager(const SPEManager &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            const SPEManager & operator= (const SPEManager &);

            /// \}

        public:
            friend class SPEInstruction;
            friend class SPEInstructionQueue;
            friend class SPEInstructionStream;

            /// Return the only instance of SPEManager.
            static SPEManager * instance();

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

            // Dispatch all SPEInstructions in the stream to one single SPE.
            void dispatch(SPEInstructionStream & instruction_stream);

            // Wait for a SPEinstructin to be finished.
            void wait(SPEInstruction & instruction);

            // Wait for all SPEInstrictions in the queue to finish.
            void wait(SPEInstructionQueue & instruction_queue);

            // Wait for all SPEInstrictions in the stream so far to finish.
            void wait(SPEInstructionStream & instruction_stream);
    };
}

#endif
