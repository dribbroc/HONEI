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
#include <libutil/spe_kernel_manager.hh>
#include <libutil/spe_manager.hh>

#include <string>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

namespace
{
    /// Compare two SPEs in terms of instruction_load
    static inline bool compare_by_load(const honei::SPE a, const honei::SPE b)
    {
        return a.kernel()->instruction_load() < b.kernel()->instruction_load();
    }

    /// Compare two SPEs in terms of last usage time
    static inline bool compare_by_last_finished(const honei::SPE a, const honei::SPE b)
    {
        return a.kernel()->last_finished() < b.kernel()->last_finished();
    }
}

namespace honei
{
    struct SPEManager::Implementation
    {
        typedef std::vector<SPE> SPEList;

        /// Our mutex.
        Mutex * const mutex;

        /// Our number of available SPEs.
        unsigned spe_count;

        /// Out list of SPEs.
        SPEList spe_list;

        /// Dispatch an SPEInstruction to an SPE.
        inline void dispatch(const SPEInstruction & instruction)
        {
            CONTEXT("When dispatching Instruction to an SPE");

            //todo
            //schauen ob ein freier kernel op unterstuetzt.
            //ansonsten einen kernel per LRU rauswerfen und neuen kernel der op beherrscht einladen.
            //wenn alle kernel beschaeftigt sind, evtl kann ein kernel op der noch am rechnen ist?
            //sonst einfach blockieren und warten bis ein kernel frei wird,
            //wenn 6 spe's rechnen koennte der user eh nicht mehr viel mehr machen

            ////Load balanced dispatching
            sort(spe_list.begin(), spe_list.end(), compare_by_load);

            LOGMESSAGE(ll_minimal, "Dispatching to SPE #" + stringify(spe_list.begin()->id()) +
                    " (load = " + stringify(spe_list.begin()->kernel()->instruction_load()) + ")");
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

        unsigned count(_imp->spe_count);
        while(count-- > 0)
        {
            SPE spe;
            _imp->spe_list.push_back(spe);
        }

        SPEKernelManager::ListIterator reference(SPEKernelManager::instance()->find("reference"));
        if (SPEKernelManager::instance()->end() == reference)
            throw InternalError("Eek!");

        for (std::vector<SPE>::iterator i(_imp->spe_list.begin()) ; i != _imp->spe_list.end() ; i++)
        {
            i->run(SPEKernel(*reference));
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

    void
    SPEManager::dispatch(std::vector<SPEInstruction>::iterator begin, std::vector<SPEInstruction>::iterator end)
    {
        Lock l(*_imp->mutex);
        for ( ; begin != end ; ++begin)
        {
            _imp->dispatch(*begin);
        }
    }

    void
    SPEManager::wait(SPEInstruction & instruction)
    {
        instruction.wait();
    }

    void
    SPEManager::wait(std::vector<SPEInstruction>::iterator begin, std::vector<SPEInstruction>::iterator end)
    {
        for ( ; begin != end ; ++begin)
        {
            (*begin).wait();
        }
    }
}
