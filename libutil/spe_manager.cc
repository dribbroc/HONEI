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

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

using namespace honei;

extern "C"
{
    extern spe_program_handle_t kernel_reference;
}

extern const Environment kernel_reference_environment;

struct SPEManager::Implementation
{
    typedef std::vector<SPE> SPEList;

    /// Our mutex.
    Mutex * const mutex;

    /// Our number of available SPEs.
    unsigned spe_count;

    /// Out list of SPEs.
    SPEList spe_list;

    /// Compare two SPE's in terms of instruction_load
    static inline bool cmp_load(const SPE a, const SPE b)
    {
        return a.kernel()->instruction_load() < b.kernel()->instruction_load();
    }

    /// Dispatch an SPEInstruction to an SPE.
    inline void dispatch(const SPEInstruction & instruction)
    {
        CONTEXT("SPEManager: When dispatching Instruction to an SPE");

        //Load balanced dispatching
        sort(spe_list.begin(), spe_list.end(), cmp_load);
        std::cout << "SPEManager: Dispatching to spe " << spe_list.begin()->id() << " with load " << spe_list.begin()->kernel()->instruction_load() << std::endl;
        instruction.enqueue_with(spe_list.begin()->kernel());
    }

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

SPEManager::SPEManager() :
    _imp(new Implementation)
{
    Lock l(*_imp->mutex);

    //unsigned count(_imp->spe_count);
    unsigned count (1);
    while(count-- > 0)
    {
        SPE spe;
        _imp->spe_list.push_back(spe);
    }

    for (std::vector<SPE>::iterator i(_imp->spe_list.begin()) ; i != _imp->spe_list.end() ; i++)
    {
        i->run(SPEKernel(kernel_reference, &kernel_reference_environment));
    }
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

