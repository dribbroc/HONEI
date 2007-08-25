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
#include <libutil/log.hh>
#include <libutil/spe_manager.hh>
#include <libutil/tags.hh>

#include <list>
#include <map>
#include <set>
#include <vector>

#include <libspe2.h>

using namespace pg512;

/**
 * \todo Implement locking for CellBackend::Implementation::SPEData
 * \fixme Allow for transfers > 16kB
 */

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

namespace
{
    inline unsigned block_size(unsigned size)
    {
        unsigned result(size), modulo(size & 0xF);

        if (modulo > 0)
            result += 16 - modulo;

        return result;
    }

    inline bool aligned(const void * address)
    {
        return 0 == (reinterpret_cast<const unsigned long>(address) & 0xF);
    }
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

    std::string print(const AllocationState & s)
    {
        switch (s)
        {
            case as_free: return "free"; break;
            case as_used: return "used"; break;
            case as_reserved: return "reserved"; break;
            default: return "bug!";
        }
    }

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

    struct SPEData;

    /// \name Convenience typedefs
    /// \{
    typedef std::set<MemoryId> IdSet;
    typedef std::map<MemoryId, CellBackend::Chunk *> ChunkMap;
    typedef std::list<AllocationInfo> AllocationList;
    typedef std::vector<SPEData *> SPEVector;
    /// \}

    /**
     * SPEData holds all data that is uniquely linked to one SPE.
     */
    struct SPEData
    {
        /// Our set of known memory ids.
        CellBackend::Implementation::IdSet known_ids;

        /// Our map of memory ids to SPE local store memory chunks.
        ChunkMap chunk_map;

        /// Our allocations.
        AllocationList allocation_list;

        /// Our SPE.
        SPE spe;

        SPEData(const SPE & s) :
            spe(s)
        {
        }
    };


    /// Our mutex.
    Mutex * const mutex;

    /// Our data.
    SPEVector spes;

    Implementation() :
        mutex(new Mutex)
    {
        for (SPEManager::Iterator i(SPEManager::instance()->begin()), i_end(SPEManager::instance()->end()) ;
                i != i_end ; ++i)
        {
            SPEData * data(new SPEData(*i));
            // \todo Get information on reserved space! Use SPEProgram, ELF information.
            // \todo Hardcoding 0x4000@0 and 0x1000@end for now.
            AllocationInfo program(0x0, 0x4000, 0x4000, as_reserved);
            AllocationInfo space(0x4000, 0x3b000, 0x3b000, as_free);
            AllocationInfo end(0x3f000, 0x1000, 0x1000, as_reserved);

            data->allocation_list.push_back(program);
            data->allocation_list.push_back(space);
            data->allocation_list.push_back(end);
            spes.push_back(data);
        }
    }

    ~Implementation()
    {
        delete mutex;

        for (SPEVector::iterator i(spes.begin()), i_end(spes.end()) ; i != i_end ; ++i)
        {
            ChunkMap & map((*i)->chunk_map);
            for (ChunkMap::iterator c(map.begin()), c_end(map.end()) ; c != c_end ; ++c)
            {
                delete c->second;
            }
            delete *i;
        }
    }

    CellBackend::Chunk * alloc(const DeviceId device, unsigned int size)
    {
        CONTEXT("When allocating space of size '" + stringify(size) + "':");
        CellBackend::Chunk * result(0);

        ASSERT(device < spes.size(), "invalid SPE id '" + stringify(device) + "'!");
        AllocationList & list(spes[device]->allocation_list);

        for (AllocationList::iterator i(list.begin()), i_end(list.end()) ; i != i_end ; ++i)
        {
            if (as_free != i->state)
                continue;

            if (size > i->block_size)
                continue;

            AllocationInfo info(i->address, size, ::block_size(size), as_used);

            i->address += info.block_size;
            i->block_size -= info.block_size;
            i->size = i->block_size;

            list.insert(i, info);

            result = new Chunk(device, info.address, info.size);

            break;
        }

        if (! result)
            throw OutOfMemoryError(tags::tv_cell, device);

        return result;
    }

    inline void swap(const MemoryId left, const MemoryId right, const DeviceId device)
    {
        CONTEXT("When swapping memory ids on SPE '" + stringify(device) + "':");
        if (left == right)
        {
            Log::instance()->message(ll_minimal, "CellBackend::swap called for identical memory ids");
            return;
        }

        ASSERT(device < spes.size(), "invalid SPE id '" + stringify(device) + "'!");
        ChunkMap & map(spes[device]->chunk_map);

        ChunkMap::iterator l(map.find(left)), r(map.find(right));
        if (map.end() == l)
            throw MemoryIdNotKnown(left, device);

        if (map.end() == r)
            throw MemoryIdNotKnown(right, device);

        Chunk * tmp(l->second);
        l->second = r->second;
        r->second = tmp;
    }

    inline void free(const DeviceId device, unsigned address, unsigned size)
    {
        CONTEXT("When freeing allocated space of size '" + stringify(size) + "' at '" + stringify(address) + "':");
        ASSERT(device < spes.size(), "invalid SPE id '" + stringify(device) + "'!");

        AllocationList & list(spes[device]->allocation_list);

        AllocationInfo info(address, size, ::block_size(size), as_used);
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
            iter->state = as_free;

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

    inline bool id_known(const MemoryId id, const DeviceId device)
    {
        IdSet & set(spes[device]->known_ids);

        return set.end() != set.find(id);
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

    if (! _imp->id_known(id, device))
        throw MemoryIdNotKnown(id, device);

    ASSERT(device < _imp->spes.size(), "invalid SPE id '" + stringify(device) + "'!");
    Implementation::ChunkMap & map(_imp->spes[device]->chunk_map);
    Implementation::ChunkMap::iterator c(map.find(id));
    ASSERT(map.end() != c, "No chunk found for memory id '" + stringify(id) + "'!");

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
    if (! aligned(address))
        throw MisalignmentError(address, 16, tags::tv_cell, device);

    Lock l(*_imp->mutex);

    ASSERT(device < _imp->spes.size(), "invalid SPE id '" + stringify(device) + "'!");
    Implementation::SPEData & spe(*_imp->spes[device]);
    spe.known_ids.insert(id);
    Implementation::ChunkMap::iterator c(spe.chunk_map.find(id)), c_end(spe.chunk_map.end());
    if (c_end != c) // We already uploaded for this memory id
    {
        if (size != c->second->size) // We uploaded a different size \todo checkfor > ?
        {
            _imp->free(c->second);
            c = spe.chunk_map.insert(spe.chunk_map.end(), std::make_pair(id, _imp->alloc(device, size)));
        }
    }
    else
    {
        c = spe.chunk_map.insert(spe.chunk_map.end(), std::make_pair(id, _imp->alloc(device, size)));
    }
    ASSERT(c_end != c, "urgh... something very wrong in CellBackend::upload()!");

    unsigned int tag(_imp->next_upload_tag()), tag_status(0), counter(5);
    do
    {
        int retval(spe_mfcio_get(spe.spe.context(), c->second->address, address, size, tag, 0, 0));
        spe_mfcio_tag_status_read(spe.spe.context(), 0, SPE_TAG_ANY, &tag_status);
        --counter;
    }
    while ((! (tag_status && (1 << tag))) && (counter > 0));
    ASSERT(counter > 0, "DMA transfer didn't work!");
}

void
CellBackend::download(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo Does not yet work for size > 16kB.
    CONTEXT("When downloading '" + stringify(size) + "' bytes of data to '" + stringify(address) +
            "' for memory id '" + stringify(id) + "' from SPE '" + stringify(device) + "':");
    if (! aligned(address))
        throw MisalignmentError(address, 16, tags::tv_cell, device);

    Lock l(*_imp->mutex);
    ASSERT(device < _imp->spes.size(), "invalid SPE id '" + stringify(device) + "'!");

    if (! _imp->id_known(id, device))
        throw MemoryIdNotKnown(id, device);

    CellBackend::Chunk * result(0);
    Implementation::SPEData & spe(*_imp->spes[device]);
    Implementation::ChunkMap::const_iterator c(spe.chunk_map.find(id)), c_end(spe.chunk_map.end());
    ASSERT(c_end != c, "No chunk found for memory id '" + stringify(id) + "'!");

    unsigned int tag(_imp->next_download_tag()), tag_status(0x3), counter(5);
    do
    {
        int retval(spe_mfcio_put(spe.spe.context(), c->second->address, address, size, tag, 0, 0));
        spe_mfcio_tag_status_read(spe.spe.context(), 0, SPE_TAG_ANY, &tag_status);
        --counter;
    }
    while ((! (tag_status && (1 << tag))) && (counter > 0));
}

void
CellBackend::swap(const MemoryId left, const MemoryId right, const DeviceId device)
{
    CONTEXT("When swapping memory ids '" + stringify(left) + "' and '" + stringify(right) + "' on " +
            "device '" + stringify(device) + "':");
    Lock l(*_imp->mutex);

    _imp->swap(left, right, device);
}

void
CellBackend::swap(const MemoryId left, const MemoryId right)
{
    CONTEXT("When swapping memory ids '" + stringify(left) + "' and '" + stringify(right) + "':");
    Lock l(*_imp->mutex);

    for (SPEManager::Iterator i(SPEManager::instance()->begin()), i_end(SPEManager::instance()->end()) ;
            i != i_end ; ++i)
    {
        bool left_known(_imp->id_known(left, i->id())), right_known(_imp->id_known(right, i->id()));
        if (left_known != right_known)
        {
            if (! left_known)
                throw MemoryIdNotKnown(left, i->id());
            else
                throw MemoryIdNotKnown(right, i->id());
        }

        if (! left_known)
            continue;

        _imp->swap(left, right, i->id());
    }
}

void
CellBackend::free(const MemoryId id, const DeviceId device)
{
    CONTEXT("When freeing memory id '" + stringify(id) + "' on SPE '" + stringify(device) + "':");
    Lock l(*_imp->mutex);

    if (! _imp->id_known(id, device))
        throw MemoryIdNotKnown(id, device);

    ASSERT(device < _imp->spes.size(), "invalid SPE id '" + stringify(device) + "'!");
    Implementation::SPEData & spe(*_imp->spes[device]);
    spe.known_ids.erase(id);
    Implementation::ChunkMap::iterator c(spe.chunk_map.find(id)), c_end(spe.chunk_map.end());
    ASSERT(c_end != c, "no chunk found for memory id '" + stringify(id) + "'!");

    _imp->free(c->second);
    spe.chunk_map.erase(c);
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
