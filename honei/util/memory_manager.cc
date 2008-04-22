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

#include <honei/util/assertion.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/memory_manager.hh>

using namespace honei;

template class InstantiationPolicy<MemoryManager, Singleton>;

MemoryAddressNotKnown::MemoryAddressNotKnown(const void * address) :
    Exception("No memory id found for address '" + stringify(address) + "'")
{
}

MemoryIdNotKnown::MemoryIdNotKnown(const MemoryId id) :
    Exception("Memory id '" + stringify(id) + "' is unknown to the memory manager")
{
}

MemoryIdNotKnown::MemoryIdNotKnown(const MemoryId id, const DeviceId device) :
    Exception("Memory id '" + stringify(id) + "' is unknown to backend device '" + stringify(device) + "'")
{
}

MemoryChunkSizeInvalid::MemoryChunkSizeInvalid(unsigned long size, const std::string & msg) :
    Exception("Chunk size '" + stringify(size) + "' is invalid: " + msg)
{
}

MemoryId
MemoryManager::associate(void * address, const std::ptrdiff_t size)
{
    CONTEXT("When associating with memory chunk of size '" + stringify(size)
            + "' at address '" + stringify(address) + "':");
    ASSERT(0 != address, "Invalid address 'NULL' for associate!");
    ASSERT(0 != size, "Invalid size '0' for associate!");

    MemoryId result(_get_free_id());
    MemoryInfo * info(new MemoryInfo(address, size, tags::tv_cpu));
    _address_map[result] = address;
    _info_map[result] = info;

    return result;
}

void
MemoryManager::upload(const MemoryId id, const DeviceId device, const tags::TagValue destination)
{
    CONTEXT("When uploading complete chunk to '" + stringify(destination)
            + "'-memory for memory id '" + stringify(id) + "':");
    ASSERT(tags::tv_cpu != destination, "Invalid destination 'CPU' for upload!");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    InfoMap::iterator i(_info_map.find(id));
    ASSERT(_info_map.end() != i, "No info map entry found for id '" + stringify(id) + "'!");

    _get_backend(destination)->upload(id, device, i->second->_address, i->second->size);

    MemoryInfo * info(new MemoryInfo(i->second->_address, i->second->size, destination));
    delete i->second;
    i->second = info;
}

void
MemoryManager::upload(const MemoryId id, const DeviceId device, const tags::TagValue destination, const std::ptrdiff_t size)
{
    CONTEXT("When uploading '" + stringify(size) + "' bytes to '" + stringify(destination)
            + "'-memory for memory id '" + stringify(id) + "':");
    ASSERT(tags::tv_cpu != destination, "Invalid destination 'CPU' for upload!");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    AddressMap::const_iterator a(_address_map.find(id));
    ASSERT(_address_map.end() != a, "No address map entry found for id '" + stringify(id) + "'!");

    _get_backend(destination)->upload(id, device, a->second, size);

    MemoryInfo * info(new MemoryInfo(a->second, size, destination));
    _info_map[id] = info;
}

void
MemoryManager::download(const MemoryId id, const DeviceId device)
{
    CONTEXT("When downloading device memory for id '" + stringify(id) + "':");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    InfoMap::const_iterator i(_info_map.find(id));
    ASSERT(_info_map.end() != i, "No info map entry found for id '" + stringify(id) + "'!");

    download(id, device, i->second->_address, i->second->size);
}

void
MemoryManager::download(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    CONTEXT("When downloading device memory for id '" + stringify(id) + "':");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    InfoMap::iterator i(_info_map.find(id));
    ASSERT(_info_map.end() != i, "No info map entry found for id '" + stringify(id) + "'!");

     _get_backend(i->second->location)->download(id, device, address, size);

    _address_map[id] = address;

    tags::TagValue t(i->second->location);
    delete i->second;
    i->second = new MemoryInfo(address, size, t);
}

void
MemoryManager::free(const MemoryId id)
{
    CONTEXT("When freeing device memory for id '" + stringify(id) + "':");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    InfoMap::iterator i(_info_map.find(id));
    ASSERT(_info_map.end() != i, "No map entry found for id '" + stringify(id) + "'!");

    _get_backend(i->second->location)->free(id);

    _address_map.erase(_address_map.find(id));
    delete i->second;
    _info_map.erase(i);
    _ids.erase(id);
}

const MemoryId
MemoryManager::get_id_by_address(const void * address) const
{
    CONTEXT("When obtaining memory id for address '" + stringify(address) + "':");

    AddressMap::const_iterator a(_address_map.begin()), a_end(_address_map.end());
    for ( ; a != a_end ; ++a)
    {
        if (address == a->second)
            break;
    }

    if (a == a_end)
        throw MemoryAddressNotKnown(address);

    return a->first;
}

MemoryInfo
MemoryManager::get_info_by_id(const MemoryId id) const
{
    CONTEXT("When obtaining information on memory id '" + stringify(id) + "':");

    if (! _id_used(id))
        throw MemoryIdNotKnown(id);

    InfoMap::const_iterator i(_info_map.find(id)), i_end(_info_map.end());
    ASSERT(i_end != i, "no information found for memory id '" + stringify(id) + "':");

    return *i->second;
}

void
MemoryManager::swap(const MemoryId left, const MemoryId right)
{
    CONTEXT("When swapping device memory information for ids '" + stringify(left) + "' and '" +
            stringify(right) + "':");

    if (! _id_used(left))
        throw MemoryIdNotKnown(left);

    if (! _id_used(right))
        throw MemoryIdNotKnown(right);

    if (_info_map[left]->size != _info_map[right]->size)
        throw InternalError("Trying to swap memory blocks of different sizes");

    if (_info_map[left]->location != _info_map[right]->location)
        throw InternalError("Trying to swap memory blocks that live on different devices");

    /// \todo This also copies transfer times. Is that correct?
    MemoryInfo * tmp_info(_info_map[left]);
    _info_map[left] = _info_map[right];
    _info_map[right] = tmp_info;

    void * tmp_address(_address_map[left]);
    _address_map[left] = _address_map[right];
    _address_map[right] = tmp_address;

    /// \todo Call swap for _all_ backends ?
    _get_backend(tmp_info->location)->swap(left, right);
}

