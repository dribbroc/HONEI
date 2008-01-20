/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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

#ifndef LIBUTIL_GUARD_SPE_TRANSFER_LIST_HH
#define LIBUTIL_GUARD_SPE_TRANSFER_LIST_HH 1

#include <honei/cell/types.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/stringify.hh>

namespace honei
{
    /**
     * TransferListSizeExceeded is thrown when an element is added although the list
     * containts max_size elements.
     *
     * \ingroup grpcell
     * \ingroup grpexceptions
     */
    struct TransferListSizeExceeded :
        public Exception
    {
        TransferListSizeExceeded() :
            Exception("Specified maximum number of list elements exceeded")
        {
        }
    };

    struct TransferListTransferSizeExceeded :
        public Exception
    {
        TransferListTransferSizeExceeded() :
            Exception("Specified maximum size of transfer list exceeded")
        {
        }
    };

    /**
     *
     * \ingroup grpcell
     */
    class SPETransferList
    {
        public:
            typedef cell::ListElement ListElement;

        private:
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

        public:
            /**
             * Constructor.
             *
             * Creates a list of transfer elements suitable for SPE-side DMA
             * list transfers.
             *
             * \param max_size Maximal number of list elements for this list.
             * \param max_transfer_size Maximal number of bytes which shall be
             *                          transfered by this list.
             */
            SPETransferList(unsigned max_size, unsigned max_transfer_size);

            /// Add an address/size pair to the list.
            ListElement * add(void * address, unsigned size, bool stall_and_notify = false);

            /// Return the effective address of our first element.
            unsigned long long effective_address() const;

            /// Return a pointer to the list's elements.
            ListElement * elements() const;

            /// Return our number of elements.
            unsigned size() const;

            /// Return our list's total transfer size.
            unsigned transfer_size() const;
    };
}

#endif
