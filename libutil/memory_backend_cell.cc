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

#include <libutil/lock.hh>
#include <libutil/memory_manager.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/mutex.hh>
#include <libutil/lock.hh>
#include <libutil/spe_manager.hh>
#include <libutil/tags.hh>

#include <list>
#include <map>
#include <set>
#include <iostream>

#include <libspe2.h>

using namespace pg512;

CellBackend::Chunk::Chunk(DeviceId i, unsigned int a, unsigned int s) :
    id(i),
    address(a),
    size(s)
{
}

CellBackend::Chunk::Chunk(const CellBackend::Chunk & other) :
    id(other.id),
    address(other.address),
    size(other.size)
{
}

CellBackend::Chunk::~Chunk()
{
}

struct CellBackend::Implementation
{
    /**
     * AllocationState describes if a memory region is either free, used by any chunk
     * or reserved.
     */
    enum AllocationState
    {
        as_free = 0,
        as_used,
        as_reserved
    };

    /**
     * AllocationInfo provides all information to uniquely identify and handle a
     * chunk of SPE local store memory.
     */
    struct AllocationInfo
    {
        unsigned address;
        unsigned size;
        unsigned block_size;
        AllocationState state;

        AllocationInfo(unsigned a, unsigned s, unsigned bs, AllocationState st) :
            address(a),
            size(s),
            block_size(bs),
            state(st)
        {
        }

        bool operator== (const AllocationInfo & other)
        {
            return (address == other.address) && (size == other.size) && (block_size == other.block_size)
                && (state == other.state);
        }
    };

    /// \name Convenience typedefs
    /// \{
    typedef std::set<MemoryId> IdSet;
    typedef std::multimap<MemoryId, CellBackend::Chunk *> ChunkMap;
    typedef std::map<DeviceId, SPE> SPEMap;
    typedef std::list<AllocationInfo> AllocationList;
    typedef std::map<DeviceId, AllocationList *> AllocationMap;
    /// \}

    /// Our set of known memory ids.
    IdSet known_ids;

    /// Our map of memory ids to SPE local store memory chunks.
    ChunkMap chunk_map;

    /// Our SPEs.
    SPEMap spe_map;

    /// Our SPEs' allocations.
    AllocationMap allocation_map;

    /// Our mutex.
    Mutex * const mutex;

    Implementation() :
        mutex(new Mutex)
    {
        for (SPEManager::Iterator i(SPEManager::instance()->begin()), i_end(SPEManager::instance()->end()) ;
                i != i_end ; ++i)
        {
            spe_map.insert(std::make_pair(i->id(), *i));
            // \todo Get information on reserved space! Use SPEProgram, ELF information.
            // \todo Hardcoding 0x4000@0 and 0x1000@end for now.
            AllocationList * list(new AllocationList);
            AllocationInfo program(0x0, 0x4000, 0x4000, as_reserved);
            AllocationInfo space(0x4000, 0x3b000, 0x3b000, as_free);
            AllocationInfo end(0x3f000, 0x1000, 0x1000, as_reserved);

            list->push_back(program);
            list->push_back(space);
            list->push_back(end);
            allocation_map[i->id()] = list;
        }
    }

    ~Implementation()
    {
        delete mutex;

        for (AllocationMap::iterator m(allocation_map.begin()), m_end(allocation_map.end()) ;
                m != m_end ; ++m)
        {
            delete m->second;
        }
    }

    CellBackend::Chunk * alloc(const DeviceId id, unsigned int size)
    {
        CONTEXT("When allocating space of size '" + stringify(size) + "':");
        CellBackend::Chunk * result(0);

        AllocationMap::iterator li(allocation_map.find(id));
        ASSERT(allocation_map.end() != li, "invalid SPE id '" + stringify(id) + "'!");
        AllocationList & list(*li->second);

        for (AllocationList::iterator i(list.begin()), i_end(list.end()) ; i != i_end ; ++i)
        {
            if (as_free != i->state)
                continue;

            if (size > i->block_size)
                continue;

            AllocationInfo info(i->address, size, ((size / 16) + 1) * 16, as_used);

            i->address += info.block_size;
            i->block_size -= info.block_size;
            i->size = i->block_size;

            list.insert(i, info);

            result = new Chunk(id, info.address, info.size);

            break;
        }

        if (! result)
            throw OutOfMemoryError(tags::tv_cell, id);

        return result;
    }

    inline void free(const DeviceId id, unsigned address, unsigned size)
    {
        CONTEXT("When freeing allocated space of size '" + stringify(size) + "' at '" + stringify(address) + "':");
        AllocationMap::iterator li(allocation_map.find(id));
        ASSERT(allocation_map.end() != li, "invalid SPE id '" + stringify(id) + "'!");
        AllocationList & list(*li->second);

        AllocationInfo info(address, size, (((size / 16) + 1) * 16), as_used);
        AllocationList::iterator iter(std::find(list.begin(), list.end(), info));
        ASSERT(list.end() != iter, "double free or corruption!");

        AllocationList::iterator prev(iter), next(iter);
        --prev; ++next;
        if ((prev->state == as_free) && (next->state == as_free))
        {
            // Triple-merge
            prev->block_size += iter->block_size + next->block_size;
            prev->size = prev->block_size;

            list.erase(iter);
            list.erase(next);
        }
        else if (prev->state == as_free)
        {
            // Merge prev and iter
            prev->block_size += iter->block_size;
            prev->size = prev->block_size;

            list.erase(iter);
        }
        else if (next->state == as_free)
        {
            // Merge iter and next
            iter->block_size += next->block_size;
            iter->size = prev->block_size;

            list.erase(next);
        }
        else
        {
            iter->state = as_free;
        }
    }

    inline void free(const CellBackend::Chunk * chunk)
    {
        CONTEXT("When freeing chunk '" + stringify(*chunk) + "':");

        free(chunk->id, chunk->address, chunk->size);

        delete chunk;
    }

    inline unsigned int next_download_tag()
    {
        static unsigned int last_download_tag(0);

        ++last_download_tag;

        return last_download_tag &= 0x1F; // Clamp tags to 0..31!
    }

    inline unsigned int next_upload_tag()
    {
        static unsigned int last_upload_tag(0);

        ++last_upload_tag;

        return last_upload_tag &= 0x1F; // Clamp tags to 0..31!
    }
};

CellBackend::CellBackend() :
    _imp(new CellBackend::Implementation)
{
}

CellBackend::~CellBackend()
{
    for (std::map<MemoryId, CellBackend::Chunk *>::iterator c(_imp->chunk_map.begin()),
            c_end(_imp->chunk_map.end()) ; c != c_end ; ++c)
    {
        free(c->second);
    }
}

CellBackend *
CellBackend::instance()
{
    static CellBackend result;

    return &result;
}

const CellBackend::Chunk *
CellBackend::get_chunk(const MemoryId id, const DeviceId device)
{
    CONTEXT("When retrieving Chunk for memory id '" + stringify(id) + "' on SPE '" + stringify(device) + ":");
    Lock l(*_imp->mutex);

    if (_imp->known_ids.end() == _imp->known_ids.find(id))
        throw MemoryIdNotKnown(id);

    Implementation::ChunkMap::const_iterator c(_imp->chunk_map.find(id)), c_end(_imp->chunk_map.upper_bound(id));
    for ( ; c != c_end ; ++c)
    {
        if (device != c->second->id)
            continue;

        break;
    }
    ASSERT(c_end != c, "No information found for id '" + stringify(id) + "'!");

    return c->second;
}

MemoryBackend *
CellBackend::backend_instance()
{
    return instance();
}

void
CellBackend::upload(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo Does not yet work for size > 16kB.
    CONTEXT("When uploading '" + stringify(size) + "' bytes of data from '" + stringify(address) +
            "' for memory id '" + stringify(id) + "' to SPE '" + stringify(device) + "':");
    Lock l(*_imp->mutex);

    _imp->known_ids.insert(id);
    Implementation::ChunkMap::iterator c(_imp->chunk_map.find(id)), c_end(_imp->chunk_map.upper_bound(id));
    for ( ; c != c_end ; ++c)
    {
        if (c->second->id != device)
            continue;

        break;
    }

    if (c_end != c) // We already uploaded for this memory id
    {
        if (size != c->second->size) // We uploaded a different size \todo checkfor > ?
        {
            _imp->free(c->second);
            c = _imp->chunk_map.insert(_imp->chunk_map.end(), std::make_pair(id, _imp->alloc(device, size)));
        }
    }
    else
    {
        c = _imp->chunk_map.insert(_imp->chunk_map.end(), std::make_pair(id, _imp->alloc(device, size)));
    }
    ASSERT(c_end != c, "urgh... something very wrong in CellBackend::upload()!");

    Implementation::SPEMap::const_iterator s(_imp->spe_map.find(c->second->id));
    ASSERT(_imp->spe_map.end() != s, "No SPE found with device id '" + stringify(c->second->id) + "'!");

    unsigned int tag(_imp->next_upload_tag()), tag_status(3), counter(5);
    do
    {
        int retval(spe_mfcio_get(s->second.context(), c->second->address, address, size, tag, 0, 0));
        spe_mfcio_tag_status_read(s->second.context(), 0, SPE_TAG_IMMEDIATE, &tag_status);
        --counter;
    }
    while ((tag_status != 0) && (counter > 0));
}

void
CellBackend::download(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo Does not yet work for size > 16kB.
    CONTEXT("When downloading '" + stringify(size) + "' bytes of data to '" + stringify(address) +
            "' for memory id '" + stringify(id) + "' from SPE '" + stringify(device) + "':");
    Lock l(*_imp->mutex);

    if (_imp->known_ids.end() == _imp->known_ids.find(id))
        throw MemoryIdNotKnown(id);

    CellBackend::Chunk * result(0);
    Implementation::ChunkMap::const_iterator c(_imp->chunk_map.find(id)), c_end(_imp->chunk_map.upper_bound(id));
    for ( ; c != c_end ; ++c)
    {
        if (device != c->second->id)
            continue;

        break;
    }
    ASSERT(c_end != c, "No chunk found for memory id '" + stringify(id) + "'!");

    Implementation::SPEMap::const_iterator s(_imp->spe_map.find(c->second->id));
    ASSERT(_imp->spe_map.end() != s, "No SPE found with device id '" + stringify(c->second->id) + "'!");

    unsigned int tag(_imp->next_download_tag()), tag_status(0), counter(5);
    do
    {
        spe_mfcio_put(s->second.context(), c->second->address, address, size, tag, 0, 0);
        spe_mfcio_tag_status_read(s->second.context(), tag, SPE_TAG_IMMEDIATE, &tag_status);
    }
    while ((tag_status != 0) && (counter > 0));
}

void
CellBackend::swap(const MemoryId left, const MemoryId right)
{
    CONTEXT("When swapping memory ids '" + stringify(left) + "' and '" + stringify(right) + "':");
    Lock l(*_imp->mutex);

    Implementation::ChunkMap::iterator i(_imp->chunk_map.find(left)), j(_imp->chunk_map.find(right));
    if (i == _imp->chunk_map.end())
        throw MemoryIdNotKnown(left);

    if (j == _imp->chunk_map.end())
        throw MemoryIdNotKnown(right);

    CellBackend::Chunk * tmp(i->second);
    i->second = j->second;
    j->second = tmp;
}

void
CellBackend::free(const MemoryId id, const DeviceId device)
{
    CONTEXT("When freeing memory id '" + stringify(id) + "':");
    Lock l(*_imp->mutex);

    if (_imp->known_ids.end() != _imp->known_ids.find(id))
        throw MemoryIdNotKnown(id);

    _imp->known_ids.erase(id);
    Implementation::ChunkMap::iterator c(_imp->chunk_map.find(id)), c_end(_imp->chunk_map.upper_bound(id));
    for ( ; c != c_end ; ++c)
    {
        if (device != c->second->id)
            continue;

        break;
    }
    ASSERT(c_end != c, "no chunk found for memory id '" + stringify(id) + "'!");
    _imp->free(c->second);
}

const CellBackend::Chunk *
CellBackend::alloc(const DeviceId device, const unsigned int size)
{
    return _imp->alloc(device, size);
}

inline void
CellBackend::free(const CellBackend::Chunk * chunk)
{
    CONTEXT("When freeing memory chunk at '" + stringify(chunk) + "':");
    Lock l(*_imp->mutex);

    _imp->free(chunk);
}

static MemoryBackendRegistrator cell_backend_registrator(tags::tv_cell, &CellBackend::backend_instance);
