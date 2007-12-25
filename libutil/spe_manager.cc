/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <libutil/configuration.hh>
#include <libutil/lock.hh>
#include <libutil/log.hh>
#include <libutil/mutex.hh>
#include <libutil/spe_instruction.hh>
#include <cell/cell.hh>
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

        /// Our last 8 used OpCodes
        cell::OpCode opcode_history[8];

        /// Our counter to the next history position
        unsigned int next_history_pos;

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
            //
            /// \todo was ist, wenn kein kernel den opcode unterstuetzt ?

            sort(spe_list.begin(), spe_list.end(), compare_by_load);

            //wait until one spe is idle
            while (spe_list.begin()->kernel()->instruction_load() != 0)
            {
                sort(spe_list.begin(), spe_list.end(), compare_by_load);
            }

            SPEList::iterator spe_it(spe_list.begin());

            //search a spe that knows our instruction and is idle
            while(spe_it != spe_list.end() && !spe_it->kernel()->knows(instruction.instruction().opcode) && spe_it->kernel()->instruction_load() == 0)
            {
                spe_it ++;
            }

            //if no spe can handle our instruction, load a new kernel in an idle spe
            if(spe_it == spe_list.end() || spe_it->kernel()->instruction_load() != 0)
            {
#ifdef DEBUG
                std::string msg("SPE map before reload: \n");
                for (int i(0) ; i < spe_count ; ++i)
                {
                    msg += "SPE #" + stringify(spe_list.at(i).id()) + " : '" + stringify(spe_list.at(i).kernel()->kernel_info().name) +"'\n";
                }
                LOGMESSAGE(ll_minimal, msg);
#endif
                //LRU
                //at least one spe should have load == 0
                sort(spe_list.begin(), spe_list.end(), compare_by_last_finished);
                spe_it = spe_list.begin();
                while (spe_it->kernel()->instruction_load() != 0)
                {
                    ++spe_it;
                }

                //Halt the spe side
                SPEInstruction inst_halt(cell::oc_halt, 0);
                inst_halt.enqueue_with(spe_it->kernel());
                inst_halt.wait();

                //find a kernel that supports the most last used ops
                SPEKernelManager::MapIterator highest(SPEKernelManager::instance()->find(instruction.instruction().opcode));
                unsigned int highest_count(0);
                for (SPEKernelManager::MapIterator k_it(SPEKernelManager::instance()->find(instruction.instruction().opcode)),
                        k_end(SPEKernelManager::instance()->upper_bound(instruction.instruction().opcode)) ;
                        k_it != k_end ; ++k_it)
                {
                    unsigned int temp_count(0);
                    for (unsigned i(0) ; i < 8 ; i++)
                    {
                        for (unsigned j(0) ; j < (*k_it).capabilities.opcode_count ; ++j)
                        {
                            if ((*k_it).capabilities.opcodes[j] == opcode_history[i]);
                            {
                                temp_count++;
                                break;
                            }
                        }
                    }
                    if (temp_count > highest_count)
                    {
                        highest = k_it;
                        highest_count = temp_count;
                    }
                }
                LOGMESSAGE(ll_minimal, "Found new kernel(" + stringify ((*highest).name) +"), supporting " +
                        stringify(highest_count) + " former used opcodes.");
                SPEKernel kernel(*highest);
                spe_it->run(kernel);
#ifdef DEBUG
                msg ="SPE map after reload: \n";
                for (int i(0) ; i < spe_count ; ++i)
                {
                    msg += "SPE #" + stringify(spe_list.at(i).id()) + " : " + stringify(spe_list.at(i).kernel()->kernel_info().name) +"\n";
                }
                LOGMESSAGE(ll_minimal, msg);
#endif
            }

            LOGMESSAGE(ll_minimal, "Dispatching to SPE #" + stringify(spe_it->id()) +
                    " (kernel = " + spe_it->kernel()->kernel_info().name + ", load = " + stringify(spe_it->kernel()->instruction_load()) + ") OpCode: " + stringify(instruction.instruction().opcode));
            instruction.enqueue_with(spe_it->kernel());

            opcode_history[next_history_pos] = instruction.instruction().opcode;
            next_history_pos = (next_history_pos + 1) % 8;
        }

        /// Constructor.
        Implementation() :
            next_history_pos(0),
            mutex(new Mutex),
            spe_count(Configuration::instance()->get_value("cell::number-of-spes",
                        spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1))) /// \todo USABLE or PHYSICAL?
        {
            std::fill(opcode_history, opcode_history + 8, cell::oc_noop);
        }

        /// Destructor.
        ~Implementation()
        {
            {
                Lock l(*mutex);
                std::vector<SPE>::iterator spe_it(spe_list.begin());
                for ( ; spe_it != spe_list.end() ; ++spe_it)
                {
                    //Halt the spe side
                    SPEInstruction inst_halt(cell::oc_halt, 0);
                    inst_halt.enqueue_with(spe_it->kernel());
                    inst_halt.wait();
                }
            }
            delete mutex;
            LOGMESSAGE(ll_minimal, "All SPE's signaled halt, destoying SPEManager.");
        }
    };

    SPEManager::SPEManager() :
        _imp(new Implementation)
    {
        CONTEXT("When creating SPEManager:");
        Lock l(*_imp->mutex);

        unsigned count(_imp->spe_count);
        while(count-- > 0)
        {
            SPE spe;
            _imp->spe_list.push_back(spe);
        }

        count = _imp->spe_count;
        unsigned long kernel_count = SPEKernelManager::instance()->size();
        unsigned long spe_per_kernel = count / kernel_count;

        std::vector<SPE>::iterator spe_it(_imp->spe_list.begin());

#if 0
        for ( ; spe_it != _imp->spe_list.end() ; ++spe_it)
        {
            spe_it->run(SPEKernel(*SPEKernelManager::instance()->find("test")));
        }
#endif
        for (SPEKernelManager::ListIterator k_it(SPEKernelManager::instance()->begin()), k_end(SPEKernelManager::instance()->end()) ;
                k_it != k_end && spe_it != _imp->spe_list.end() ; ++k_it)
        {
            for (unsigned long i(0) ; i < spe_per_kernel && spe_it != _imp->spe_list.end() ; ++i)
            {
                spe_it->run(SPEKernel(*k_it));
                spe_it ++;
            }
        }
        for (SPEKernelManager::ListIterator k_it(SPEKernelManager::instance()->begin()), k_end(SPEKernelManager::instance()->end()) ;
                k_it != k_end && spe_it != _imp->spe_list.end(); ++k_it, ++spe_it)
        {
                spe_it->run(SPEKernel(*k_it));
        }
#ifdef DEBUG
        std::string msg("SPE map after constructor: \n");
        for (int i(0) ; i < count ; ++i)
        {
            msg += "SPE #" + stringify(_imp->spe_list.at(i).id()) + " : '" + stringify(_imp->spe_list.at(i).kernel()->kernel_info().name) +"'\n";
        }
        LOGMESSAGE(ll_minimal, msg);
#endif
    }

    SPEManager::~SPEManager()
    {
        CONTEXT("When destroying SPEManager:");

        delete _imp;
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
            begin->wait();
        }
    }
}
