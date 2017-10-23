/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008, 2009, 2011 Sven Mallach <mallach@honei.org>
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
#ifndef HONEI_GUARD_UTIL_TICKET_HH
#define HONEI_GUARD_UTIL_TICKET_HH 1

#include <honei/util/assertion.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/tags.hh>

#include <memory>

namespace honei
{
    struct TicketBaseImpl
    {
        Mutex mutex;

        ConditionVariable completion;

        bool completed;

        TicketBaseImpl() :
            completed(false)
        {
        }

        virtual ~TicketBaseImpl()
        {
        }

        void mark()
        {
            CONTEXT("When marking a ticket for completion:");
            Lock l(mutex);

            ASSERT(! completed, "ticket marked more than once!");

            completed = true;
            completion.signal();
        }

        void wait()
        {
            CONTEXT("When waiting for ticket completion:");
            Lock l(mutex);

            while (! completed)
            {
                completion.wait(mutex);
            }
        }
    };


    /**
     * TicketBase is the base class for all types
     * of Tickets implemented for HONEIs backends.
     * Never use directly!
     */

    // Forward declaration.
    class TicketVector;

    class TicketBase
    {
        private:

        protected:

            std::shared_ptr<TicketBaseImpl> _base;

            /// Constructor.
            TicketBase();

        public:

            /// \name Basic Operations
            /// \{

            /// Destructor.
            virtual ~TicketBase();
            /// \}

            /// Mark ticket as completed.
            /// Do not overload!
            void mark() const;

            /// Wait for ticket completion.
            /// Do not overload!
            void wait() const;

            /// \name Friends of TicketBase
            /// \{
            friend class TicketVector;
            /// \}
    };

    /**
     * Ticket is used by all asynchronous function calls to relay/query
     * information on a function's completion status.
     */

    template <typename Tag_> class Ticket;

    template <> class Ticket<tags::CPU> :
        public TicketBase,
        public PrivateImplementationPattern<Ticket<tags::CPU>, Shared>
    {
        /// \name Basic Operations
        /// \{

        /// Constructor.
        Ticket();

        virtual ~Ticket();

        /// \}

        // Derive mark and wait from TicketBase
    };

    /**
     * TicketVector is used when multiple tickets are needed.
     */

    class TicketVector :
        public InstantiationPolicy<TicketVector, NonCopyable>,
        public PrivateImplementationPattern<TicketVector, Single>
    {
        public:
            /// \name Basic Operations
            /// \{

            /// Constructor.
            TicketVector();

            /// Destructor.
            ~TicketVector();

            /// \}

            /// Push a ticket to the back of the list.
            void push_back(TicketBase ticket);

            /// Wait for ticket completion.
            void wait() const;
    };
}

#endif
