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

#include <libutil/assertion.hh>
#include <libutil/condition_variable.hh>
#include <libutil/exception.hh>
#include <libutil/lock.hh>
#include <libutil/log.hh>
#include <libutil/mutex.hh>
#include <libutil/spe_event.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_kernel.hh>
#include <libutil/spe_manager.hh>
#include <libutil/sync_point.hh>

#include <limits>
#include <sys/time.h>
#include <fstream>

namespace honei
{
    using namespace cell;

    struct SPEKernel::Implementation
    {
        /// Our SPE program.
        spe_program_handle_t handle;

        /// \name Shared data
        /// \{

        /// Our program's environment.
        Environment * environment;

        /// Our queue of instructions.
        Instruction * instructions;

        /// Our index on the currently executed instruction.
        unsigned spe_instruction_index;

        /// Our index on the next free instruction.
        unsigned next_free_index;

        /// Our counter for enqueued instructions.
        unsigned enqueued_counter;

        /// Our counter for finished instructions.
        unsigned finished_counter;

        /// Were there instructions enqueued before loading the kernel in the spe?
        bool initial_instructions;

        /// The timepoint we finished our last instruction
        timeval last_finished;

        /// \}

        /// \name Synchronization data and PPE-side kernel thread
        /// \{

        /// Our mutex.
        Mutex * const mutex;

        /// Our SPE-part-has-been-loaded synchronisation point.
        SyncPoint * kernel_loaded;

        /// Our instruction-has-finished condition variable.
        ConditionVariable * instruction_finished;

        /// Our PPE-side kernel thread.
        pthread_t * const thread;

        /// Our PPE-side kernel thread's attributes.
        pthread_attr_t * const attr;

        /// \}

        /// Our SPE.
        SPE * spe;

        static void * kernel_thread(void * argument)
        {
            LOGMESSAGE(ll_minimal, "SPEKernel: Kernel thread started (anonymous)");

            Implementation * imp(static_cast<Implementation *>(argument));

            SPE * spe(0);
            {
                Lock l(*imp->mutex);

                imp->kernel_loaded->signal_and_wait(imp->mutex);
                spe = imp->spe;
            }

            LOGMESSAGE(ll_minimal, "SPEKernel: Revceived have-been-loaded signal from SPE #" +
                    stringify(spe->id()));
            CONTEXT("When handling SPE events for SPE #" + stringify(spe->id()));

            SPEEvent event(*imp->spe, SPE_EVENT_OUT_INTR_MBOX | SPE_EVENT_SPE_STOPPED);
            Instruction current_instruction;

            do
            {
                LOGMESSAGE(ll_minimal, "SPEKernel: Waiting for event");

                spe_event_unit_t * e(event.wait(-1));

                LOGMESSAGE(ll_minimal, "SPEKernel: Event occured");
                ASSERT(! (e->events & ~(SPE_EVENT_OUT_INTR_MBOX | SPE_EVENT_SPE_STOPPED)),
                        "unexpected event happened!");

                if (e->events & SPE_EVENT_OUT_INTR_MBOX)
                {
                    LOGMESSAGE(ll_minimal, "SPEKernel: Mail event pending");
                    unsigned mail(0);

                    int retval(spe_out_intr_mbox_read(spe->context(), &mail, 1, SPE_MBOX_ALL_BLOCKING));
                    ASSERT(1 == retval, "weird return value in spe_*_mbox_read!");

                    do
                    {
                        switch (mail)
                        {
                            case km_result_dword:
                                LOGMESSAGE(ll_minimal, "SPEKernel: km_result_dword received");
                                {
                                    Lock ll(*imp->mutex);

                                    spe_out_mbox_read(spe->context(),
                                            reinterpret_cast<unsigned int *>(imp->instructions[imp->spe_instruction_index].a.ea), 1);
                                }
                                continue;

                            case km_result_qword:
                                LOGMESSAGE(ll_minimal, "SPEKernel: km_result_qword received");
                                {
                                    Lock ll(*imp->mutex);

                                    spe_out_mbox_read(imp->spe->context(),
                                            reinterpret_cast<unsigned int *>(imp->instructions[imp->spe_instruction_index].a.ea), 2);
                                }
                                continue;

                            case km_instruction_finished:
                                LOGMESSAGE(ll_minimal, "SPEKernel: km_instruction_finished received");
                                {
                                    Lock ll(*imp->mutex);

                                    imp->instructions[imp->spe_instruction_index].opcode = oc_noop;
                                    ++imp->finished_counter;
                                    ++imp->spe_instruction_index;
                                    imp->spe_instruction_index %= 8; /// \todo remove hardcoded numbers
                                    current_instruction = imp->instructions[imp->spe_instruction_index];
                                    imp->instruction_finished->broadcast();
                                    gettimeofday(&(imp->last_finished), 0);
                                }
                                continue;

                            case km_unknown_opcode:
                                throw InternalError("SPEKernel: Kernel reported invalid opcode");
                                continue;

                            case km_debug_put:
                                {
                                    Lock ll(*imp->mutex);

                                    union
                                    {
                                        EffectiveAddress ea;
                                        unsigned u[2];
                                    } ea = { 0 };
                                    LocalStoreAddress lsa(0);
                                    unsigned size(0);

                                    spe_out_mbox_read(spe->context(), ea.u, 2);
                                    spe_out_mbox_read(spe->context(), &lsa, 1);
                                    spe_out_mbox_read(spe->context(), &size, 1);
                                    LOGMESSAGE(ll_minimal, "SPEKernel: PUT transfer started, ea = " +
                                            stringify(ea.ea) + ", lsa = " + stringify(lsa) + " size = " +
                                            stringify(size));
                                }
                                continue;

                            case km_debug_get:
                                {
                                    Lock ll(*imp->mutex);

                                    union
                                    {
                                        EffectiveAddress ea;
                                        unsigned u[2];
                                    } ea = { 0 };
                                    LocalStoreAddress lsa(0);
                                    unsigned size(0);

                                    spe_out_mbox_read(spe->context(), ea.u, 2);
                                    spe_out_mbox_read(spe->context(), &lsa, 1);
                                    spe_out_mbox_read(spe->context(), &size, 1);
                                    LOGMESSAGE(ll_minimal, "SPEKernel: GET transfer started, ea = " +
                                            stringify(ea.ea) + ", lsa = " + stringify(lsa) + " size = " +
                                            stringify(size));
                                }
                                continue;

                            case km_debug_enter:
                                LOGMESSAGE(ll_minimal, "SPEKernel: SPE entering instruction");
                                continue;

                            case km_debug_leave:
                                LOGMESSAGE(ll_minimal, "SPEKernel: SPE leaving instruction");
                                continue;

                            case km_debug_acquire_block:
                                {
                                    Lock ll(*imp->mutex);

                                    LocalStoreAddress lsa(0);
                                    spe_out_mbox_read(spe->context(), &lsa, 1);
                                    LOGMESSAGE(ll_minimal, "SPEKernel: Acquired block at " +
                                            stringify(lsa));
                                }
                                continue;

                            case km_debug_release_block:
                                {
                                    Lock ll(*imp->mutex);

                                    LocalStoreAddress lsa(0);
                                    spe_out_mbox_read(spe->context(), &lsa, 1);
                                    LOGMESSAGE(ll_minimal, "SPEKernel: Released block at " +
                                            stringify(lsa));
                                }
                                continue;

                            case km_debug_dump:
                                static unsigned counter(0);
                                {
                                    Lock ll(*imp->mutex);

                                    std::string name("spu-dump-" + stringify(counter));
                                    std::fstream file(name.c_str(), std::ios::out);
                                    char * area(static_cast<char *>(spe_ls_area_get(spe->context())));
                                    for (char * c(area), * c_end(area + 256 * 1024) ; c != c_end ; ++c)
                                    {
                                        file << *c;
                                    }
                                    LOGMESSAGE(ll_minimal, "SPEKernel: Dumped LS content to file '" +
                                            name + "'");
                                    ++counter;
                                    spe->signal();
                                }
                                continue;
                        }

                        throw InternalError("Unexpected mail received in interrupt mailbox!");
                    } while (false);
                }
                else if (e->events & SPE_EVENT_SPE_STOPPED)
                {
                    LOGMESSAGE(ll_minimal, "SPEKernel: SPE stopped and signaled");
                    Lock ll(*imp->mutex);

                    imp->instruction_finished->broadcast();
                    break;
                }
                else
                {
                    LOGMESSAGE(ll_minimal, "SPEKernel: Neither mail event nor stop event");
                }
            } while (true);

            pthread_exit(0);
        }

        Implementation(const spe_program_handle_t & h, const Environment * e) :
            handle(h),
            environment(0),
            instructions(0),
            spe_instruction_index(0),
            next_free_index(0),
            enqueued_counter(0),
            finished_counter(0),
            initial_instructions(false),
            mutex(new Mutex),
            kernel_loaded(new SyncPoint(2)),
            instruction_finished(new ConditionVariable),
            thread(new pthread_t),
            attr(new pthread_attr_t),
            spe(0)
        {
            int retval(0);

            if (0 != posix_memalign(reinterpret_cast<void **>(&instructions), 128, 8 * sizeof(Instruction))) /// \todo remove hardcoded numbers
                throw std::bad_alloc(); /// \todo exception hierarchy

            if (0 != posix_memalign(reinterpret_cast<void **>(&environment), 16, sizeof(Environment)))
                throw std::bad_alloc(); /// \todo exception hierarchy

            *environment = *e;

            if (0 != (retval = pthread_attr_init(attr)))
                throw PThreadError("pthread_attr_init", retval);

            if (0 != (retval = pthread_attr_setdetachstate(attr, PTHREAD_CREATE_JOINABLE)))
                throw PThreadError("pthread_attr_setdetachstate", retval);

            if (0 != (retval = pthread_create(thread, 0, &kernel_thread, this)))
                throw PThreadError("pthread_create failed", retval);

            Instruction noop = { oc_noop };
            std::fill(instructions, instructions + 8, noop);

            gettimeofday(&(last_finished), 0);
        }

        ~Implementation()
        {
            int retval(0);

            if (0 != (retval = pthread_join(*thread, 0)))
                throw PThreadError("pthread_join", retval);

            free(instructions);
            free(environment);

            if (0 != (retval = pthread_attr_destroy(attr)))
                throw PThreadError("pthread_attr_destroy", retval);

            delete attr;
            delete thread;
            delete instruction_finished;
            delete kernel_loaded;
            delete mutex;
        }

    };

    SPEKernel::SPEKernel(const spe_program_handle_t & handle, const Environment * environment) :
        _imp(new Implementation(handle, environment))
    {
    }

    unsigned
    SPEKernel::enqueue(const SPEInstruction & instruction)
    {
        CONTEXT("When enqueueing instruction:");
        unsigned result;
        Lock l(*_imp->mutex);

        while (_imp->spe_instruction_index == (_imp->next_free_index + 1) % 8) /// \todo remove hardcoded numbers
        {
            _imp->instruction_finished->wait(*_imp->mutex);
        }

        _imp->instructions[_imp->next_free_index] = instruction.instruction();
        result = _imp->enqueued_counter;
        ++_imp->next_free_index;
        ++_imp->enqueued_counter;
        _imp->next_free_index %= 8; /// \todo remove hardcoded numbers

        if (_imp->spe)
        {
            _imp->spe->signal();
        }
        else
        {
            _imp->initial_instructions = true;
        }

        return result;
    }


    void
    SPEKernel::wait(unsigned instruction_index) const
    {
        CONTEXT("When waiting until instruction '" + stringify(instruction_index) + "' has finished:");
        Lock l(*_imp->mutex);

        while (instruction_index >= _imp->finished_counter)
        {
            _imp->instruction_finished->wait(*_imp->mutex);
        }
    }

    void
    SPEKernel::load(const SPE & spe) const
    {
        CONTEXT("When loading kernel into SPE #" + stringify(spe.id()) + ":");
        Lock l(*_imp->mutex);
        if (0 != spe_program_load(spe.context(), &_imp->handle))
        {
            throw SPEError("spe_program_load", errno);
        }

        _imp->spe = new SPE(spe);
        _imp->kernel_loaded->signal_and_wait(_imp->mutex);

        LOGMESSAGE(ll_minimal, "SPEKernel: Loading kernel into SPE");

        if (_imp->initial_instructions)
        {
            _imp->spe->signal();
            _imp->initial_instructions = false;
        }
    }

    unsigned
    SPEKernel::instruction_load() const
    {
        Lock l (*_imp->mutex);
        return (_imp->next_free_index - _imp->spe_instruction_index) % 8; /// \todo remove hardcoded numbers
    }

    timeval
    SPEKernel::last_finished() const
    {
        Lock l (*_imp->mutex);
        if ((_imp->next_free_index - _imp->spe_instruction_index) % 8 != 0)
        {
            timeval max;
            max.tv_sec = LONG_MAX;
            max.tv_usec = LONG_MAX;
            return max;
        }
        return _imp->last_finished;
    }

    void * SPEKernel::argument() const
    {
        Lock l(*_imp->mutex);

        return _imp->instructions;
    }

    void * SPEKernel::environment() const
    {
        Lock l(*_imp->mutex);

        return _imp->environment;
    }
}
