/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/ppe/spe_instruction.hh>
#include <honei/backends/cell/ppe/spe_kernel.hh>
#include <honei/backends/cell/ppe/spe_kernel_manager.hh>
#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/util/configuration.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <string>
#include <tr1/functional>

#include <libspe2.h>

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
    template class InstantiationPolicy<SPEManager, Singleton>;

    template <> struct Implementation<SPEManager>
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

        inline SPEList::iterator prepare_spe(const SPEInstruction & instruction)
        {
            CONTEXT("When preparing a SPE for instruction dispatching:");

            //check if any kernels knows the new opcode
            if (SPEKernelManager::instance()->find(instruction.instruction().opcode)
                    == SPEKernelManager::instance()->upper_bound(instruction.instruction().opcode))
                    throw InternalError("SPEManager: No kernel supports opcode: " + stringify(instruction.instruction().opcode) + ".");

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

                // Halt the spe side.
                SPEInstruction inst_halt(cell::oc_halt, 0);
                spe_it->kernel()->enqueue(inst_halt);
                inst_halt.wait();

                // Find a kernel that supports the most recently used operations.
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
                            if ((*k_it).capabilities.opcodes[j] == opcode_history[i])
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
            return spe_it;
        }
        /// Dispatch a SPEInstruction to a SPE.
        inline void dispatch(const SPEInstruction & instruction)
        {
            CONTEXT("When dispatching Instruction to a SPE:");
            SPEList::iterator spe_it(prepare_spe(instruction));
            LOGMESSAGE(ll_minimal, "Dispatching Instruction to SPE #" + stringify(spe_it->id()) +
                    " (kernel = " + spe_it->kernel()->kernel_info().name + ", load = " + stringify(spe_it->kernel()->instruction_load()) + ") OpCode: " + stringify(instruction.instruction().opcode));
            spe_it->kernel()->enqueue(instruction);

            opcode_history[next_history_pos] = instruction.instruction().opcode;
            next_history_pos = (next_history_pos + 1) % 8;
        }

        /// Dispatch a SPEInstructionQueue to a SPE
        inline void dispatch(SPEInstructionQueue & instruction_queue)
        {
            CONTEXT("When dispatching InstructionQueue to a SPE:");
            if (instruction_queue.size() > 0)
            {
                SPEList::iterator spe_it = prepare_spe(*instruction_queue.begin());
                LOGMESSAGE(ll_minimal, "Dispatching InstructionQueue to SPE #" + stringify(spe_it->id()) +
                        " (kernel = " + spe_it->kernel()->kernel_info().name + ", load = " + stringify(spe_it->kernel()->instruction_load()) + ") OpCode: " + stringify((*(instruction_queue.begin())).instruction().opcode));

                for (std::list<SPEInstruction>::iterator i(instruction_queue.begin()), i_end(instruction_queue.end()) ; i != i_end ; ++i)
                {
                    spe_it->kernel()->enqueue_queue_element(*i);
                }

                opcode_history[next_history_pos] = (*instruction_queue.begin()).instruction().opcode;
                next_history_pos = (next_history_pos + 1) % 8;
                spe_it->kernel()->run_queue();
            }
        }

        /// Constructor.
        Implementation() :
            next_history_pos(0),
            mutex(new Mutex),
            spe_count(std::min(Configuration::instance()->get_value("cell::number-of-spes",
                        spe_cpu_info_get(SPE_COUNT_USABLE_SPES, 0)), spe_cpu_info_get(SPE_COUNT_USABLE_SPES, 0))) /// \todo USABLE or PHYSICAL?
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
                    spe_it->kernel()->enqueue(inst_halt);
                    inst_halt.wait();
                }
            }
            delete mutex;
            LOGMESSAGE(ll_minimal, "All SPE's signaled halt, destoying SPEManager.");
        }
    };

    SPEManager::SPEManager() :
        PrivateImplementationPattern<SPEManager, Single>(new Implementation<SPEManager>)
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
    SPEManager::dispatch(SPEInstructionQueue & instruction_queue)
    {
        Lock l(*_imp->mutex);

        _imp->dispatch(instruction_queue);
    }

    void
    SPEManager::wait(SPEInstruction & instruction)
    {
        instruction.wait();
    }

    void
    SPEManager::wait(SPEInstructionQueue & instruction_queue)
    {
        instruction_queue.wait();
    }

    unsigned long
    SPEManager::spe_count()
    {
        return _imp->spe_count;
    }
}
