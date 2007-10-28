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
    typedef std::vector<SPE> SPEList;

    /// Our mutex.
    Mutex * const mutex;

    /// Our number of available SPEs.
    unsigned spe_count;

    /// Out list of SPEs.
    SPEList spe_list;

    /// Iterator pointing to next work-free SPE
    std::vector<SPE>::iterator next_spe;

    /// Dispatch an SPEInstruction to an SPE.
    inline void dispatch(const SPEInstruction & instruction)
    {
        CONTEXT("SPEManager: When dispatching Instruction to an SPE");

        /// \todo Find matching kernel, enqueue with that kernel
        //Round Robin dispatching
        std::cout << "SPEManager: Dispatching to spe "<< next_spe - spe_list.begin()<<" with load "<<next_spe->kernel()->instruction_load()<<std::endl;
        instruction.enqueue_with(next_spe->kernel());
        next_spe ++;
        if (next_spe == spe_list.end())
            next_spe = spe_list.begin();
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

static Environment environment __attribute__((aligned(16))) = { 0x4000, 0x3F000 };

SPEManager::SPEManager() :
    _imp(new Implementation)
{
    Lock l(*_imp->mutex);

    //unsigned count(_imp->spe_count);
    unsigned count (2);
    while(count-- > 0)
    {
        SPE spe;
        _imp->spe_list.push_back(spe);
        std::cout.flush();
    }
    for (std::vector<SPE>::iterator i(_imp->spe_list.begin()) ; i != _imp->spe_list.end() ; i++)
    {
        i->run(SPEKernel(kernel_reference, &environment));
    }
    _imp->next_spe = _imp->spe_list.begin();
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

