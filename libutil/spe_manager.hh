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

#ifndef LIBUTIL_GUARD_MEMORY_MANAGER_HH
#define LIBUTIL_GUARD_MEMORY_MANAGER_HH 1

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/memory_backend.hh>
#include <libutil/stringify.hh>

#include <libwrapiter/libwrapiter_forward_iterator.hh>

#include <string>
#include <tr1/memory>

#include <libspe2.h>

namespace pg512
{
    /**
     * SPEError is thrown by SPEManager and related classes whenever an error
     * occurs in interfacing Libspe2.
     *
     * \ingroup grpexceptions
     * \ingroup grpcell
     */
    struct SPEError :
        public ExternalError
    {
        /**
         * Constructor.
         *
         * \param msg The error message.
         * \param reason The reason for the error message.
         */
        SPEError(const std::string & msg, const std::string & reason);
    };

    class SPETask
    {
    };

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

            /// Constructor.
            SPE();

        public:
            friend class SPEManager;

            /// Copy-constructor.
            SPE(const SPE & other);

            /// Destructor.
            ~SPE();

            /// \name Iteration over our queued threads.
            /// \{

            typedef libwrapiter::ForwardIterator<SPE, SPETask> Iterator;

            Iterator begin() const;
            Iterator end() const;

            /// \}

            /// Return our device id.
            DeviceId id();
    };

    /**
     * SPEManager handles all available SPEs and dispatches SPEThread instances
     * them.
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

            /// \name Basic Operations. 
            /// \[

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
            /// Return the only instance of SPEManager.
            static SPEManager * instance();

            /// \name Iterate over our SPEs
            /// \{

            typedef libwrapiter::ForwardIterator<SPEManager, SPE> Iterator;

            Iterator begin() const;
            Iterator end() const;

            /// \}

            /// \todo dispatch 
    };
}

#endif
