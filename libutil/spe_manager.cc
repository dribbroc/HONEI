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
#include <libutil/log.hh>
#include <libutil/mutex.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_kernel.hh>
#include <libutil/spe_manager.hh>
#include <libutil/worker.hh>

#include <list>
#include <string>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

using namespace honei;

extern "C"
{
    extern spe_program_handle_t kernel_reference;
}

struct SPEManager::Implementation
{
    typedef std::list<SPE> SPEList;

    /// Our mutex.
    Mutex * const mutex;

    /// Our number of available SPEs.
    unsigned spe_count;

    /// Out list of SPEs.
    SPEList spe_list;

    /// Dispatch an SPEInstruction to an SPEs.
    inline void dispatch(const SPEInstruction & instruction)
    {
        /// \todo Find matching kernel, enqueue with that kernel
        instruction.enqueue_with(spe_list.begin()->kernel());
    }

    static Environment environment __attribute__((aligned(16)));

    /// Constructor.
    Implementation() :
        mutex(new Mutex),
        spe_count(spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1)) /// \todo USABLE or PHYSICAL?
    {
    }

    /// Destructor.
    ~Implementation()
    {
        delete mutex;
    }
};

Environment
SPEManager::Implementation::environment = { 0x4000, 0x3F000 };

SPEManager::SPEManager() :
    _imp(new Implementation)
{
    Lock l(*_imp->mutex);

    unsigned count(_imp->spe_count);
    while(count-- > 0)
    {
        SPE spe;
//        spe.run(SPEKernel(kernel_reference, &Implementation::environment));
        _imp->spe_list.push_front(spe);
    }
    _imp->spe_list.begin()->run(SPEKernel(kernel_reference, &Implementation::environment));
}

SPEManager::~SPEManager()
{
}

SPEManager *
SPEManager::instance()
{
    static SPEManager result;

    return &result;
}

SPEManager::Iterator
SPEManager::begin() const
{
    Lock l(*_imp->mutex);

    return Iterator(_imp->spe_list.begin());
}

SPEManager::Iterator
SPEManager::end() const
{
    Lock l(*_imp->mutex);

    return Iterator(_imp->spe_list.end());
}

unsigned int
SPEManager::spe_count() const
{
    return _imp->spe_count;
}

void
SPEManager::dispatch(const SPEInstruction & instruction)
{
    Lock l(*_imp->mutex);

    _imp->dispatch(instruction);
}

