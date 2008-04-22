/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/spe_transfer_list.hh>

#include <tr1/memory>

namespace honei
{
    template <> struct Implementation<SPETransferList>
    {
        typedef cell::ListElement ListElement;

        /// Our elements.
        ListElement * elements;

        /// Our maximal transfer size per element. (must be < 16 KB)
        unsigned max_transfer_size;

        /// Our maximal number of elements. (must be < 2048)
        unsigned max_size;

        /// Our actual transfer size.
        unsigned transfer_size;

        /// Our actual number of elements.
        unsigned size;

        /// Effective address (64 bit) of the last added element. They all share the upper 32 bit.
        unsigned long long effective_address;

        Implementation(unsigned s, unsigned ts) :
            effective_address(0),
            elements(new ListElement[s]),
            max_transfer_size(ts),
            max_size(s),
            transfer_size(0),
            size(0)
        {
        }

        ~Implementation()
        {
            delete[] elements;
        }
    };
}

using namespace honei;

SPETransferList::SPETransferList(unsigned max_size, unsigned max_transfer_size) :
    PrivateImplementationPattern<SPETransferList, Shared>(new Implementation<SPETransferList>(max_size, max_transfer_size))
{
    ASSERT(max_size <= 2048, "Specified number of maximum ListElements exceeds 2048");
    ASSERT(max_transfer_size <= 16384, "Specified maximum TransferList size exceeds 16 KB per ListElement");
}

SPETransferList::~SPETransferList()
{
}

SPETransferList::ListElement *
SPETransferList::add(void * address, unsigned size, bool stall_and_notify)
{
    if (size > _imp->max_transfer_size)
        throw TransferListTransferSizeExceeded();

    unsigned long long element_ea(reinterpret_cast<unsigned long long>(address));

    if ( (((element_ea >> 32) == (_imp->effective_address >> 32)) || (_imp->effective_address == 0))
            && (_imp->size + 1 <= _imp->max_size) )
    {
        _imp->effective_address = element_ea;

        ListElement * result(_imp->elements + _imp->size);
        result->stall_and_notify = stall_and_notify;
        result->reserved = 0;
        result->transfer_size = size;

        result->effective_address_low = element_ea & 0xFFFFFFFF;
        ++_imp->size;
        _imp->transfer_size += size;

        return result;
    }
    else
    {
        return 0;
    }
}

unsigned long long
SPETransferList::effective_address() const
{
    return _imp->effective_address;
}

SPETransferList::ListElement *
SPETransferList::elements() const
{
    return _imp->elements;
}

unsigned
SPETransferList::size() const
{
    return _imp->size;
}

unsigned
SPETransferList::transfer_size() const
{
    return _imp->transfer_size;
}
