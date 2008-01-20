/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Dirk Ribbrock <d_ribbrock@web.de>
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

#include <honei/libutil/lock.hh>
#include <honei/libutil/memory_manager.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/mutex.hh>
#include <honei/libutil/condition_variable.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/log.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/tags.hh>

#include <list>
#include <map>
#include <set>
#include <vector>

#include <honei/libspe2.h>

using namespace honei;



struct CellBackend::AccessCounter
{
    unsigned int _write_count;
    unsigned int _read_count;

    AccessCounter()
    {
        _write_count = 0;
        _read_count = 0;
    }
};

/**
 * \todo Implement locking for CellBackend::Implementation::SPEData
 * \fixme Allow for transfers > 16kB
 */
struct CellBackend::Implementation
{
    /// Our mutex.
    Mutex * const mutex;

    /// Our counters.
    SPECounters spe_counters;

    /// Our one access-has-finished condition variable.
    ConditionVariable * access_finished;

    Implementation() :
        mutex(new Mutex),
        spe_counters(SPEManager::instance()->spe_count()),
        access_finished(new ConditionVariable())
    {
    }

    ~Implementation()
    {
        delete mutex;
        delete access_finished;
    }

    inline void swap(const MemoryId left, const MemoryId right, const DeviceId device)
    {
        CONTEXT("When swapping memory ids on SPE '" + stringify(device) + "':");
    }

    inline void free(const DeviceId device, unsigned address, unsigned size)
    {
        CONTEXT("When freeing allocated space of size '" + stringify(size) + "' at '" + stringify(address) + "':");
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

MemoryBackend *
CellBackend::backend_instance()
{
    return instance();
}

void
CellBackend::upload(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    CONTEXT("When uploading '" + stringify(size) + "' bytes of data from '" + stringify(address) +
            "' for memory id '" + stringify(id) + "' to SPE '" + stringify(device) + "':");
}

void
CellBackend::read_enable(const DeviceId device_id)
{
    CONTEXT("When enabling reading for SPE '" + stringify(device_id) + "':");
    Lock l(*_imp->mutex);
    while (_spe_counters.at(device_id)._write_count > 0)
    {
        _imp->access_finished->wait(*_imp->mutex);
    }
    _spe_counters.at(device_id)._read_count++;
}

void
CellBackend::read_disable(const DeviceId device_id)
{
    CONTEXT("When disabling reading for SPE '" + stringify(device_id) + "':");
    Lock l(*_imp->mutex);

    if( _spe_counters.at(device_id)._read_count == 0)
    {
        //throw MemoryBackendError(tags::tv_cell, device_id, "Disabling reading impossible!");
    }
    else
    {
        _spe_counters.at(device_id)._read_count--;
        _imp->access_finished->broadcast();
    }
}


void
CellBackend::write_enable(const DeviceId device_id)
{
    CONTEXT("When enabling writing for SPE '" + stringify(device_id) + "':");
    Lock l(*_imp->mutex);

    while (_spe_counters.at(device_id)._read_count > 0 && _spe_counters.at(device_id)._read_count > 0)
    {
        _imp->access_finished->wait(*_imp->mutex);
    }
    _spe_counters.at(device_id)._write_count++;
}

void
CellBackend::write_disable(const DeviceId device_id)
{
    CONTEXT("When enabling writing for SPE '" + stringify(device_id) + "':");
    Lock l(*_imp->mutex);

    if( _spe_counters.at(device_id)._write_count == 0)
    {
        //throw MemoryBackendError("Disabling writing impossible!");
    }
    else
    {
        _spe_counters.at(device_id)._write_count--;
        _imp->access_finished->broadcast();
    }
}


void
CellBackend::download(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size)
{
    /// \todo Does not yet work for size > 16kB.
    CONTEXT("When downloading '" + stringify(size) + "' bytes of data to '" + stringify(address) +
            "' for memory id '" + stringify(id) + "' from SPE '" + stringify(device) + "':");
}

void
CellBackend::swap(const MemoryId left, const MemoryId right, const DeviceId device)
{
    CONTEXT("When swapping memory ids '" + stringify(left) + "' and '" + stringify(right) + "' on " +
            "device '" + stringify(device) + "':");
}

void
CellBackend::swap(const MemoryId left, const MemoryId right)
{
    CONTEXT("When swapping memory ids '" + stringify(left) + "' and '" + stringify(right) + "':");
}

void
CellBackend::free(const MemoryId id, const DeviceId device)
{
    CONTEXT("When freeing memory id '" + stringify(id) + "' on SPE '" + stringify(device) + "':");
}

static MemoryBackendRegistrator cell_backend_registrator(tags::tv_cell, &CellBackend::backend_instance);
