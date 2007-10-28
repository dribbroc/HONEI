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
#include <libutil/mutex.hh>
#include <libutil/spe_event.hh>
#include <libutil/spe_kernel.hh>
#include <libutil/spe_manager.hh>
#include <libutil/sync_point.hh>

#include <iostream>
namespace honei
{
    struct SPEKernel::Implementation
    {
        /// Our SPE program.
        spe_program_handle_t handle;

        /// \name Shared data
        /// \{

        /// Our program's environment.
        Environment * const environment;

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
            std::cout << "Anonymous kernel_thread startet" << std::endl;
            Implementation * imp(static_cast<Implementation *>(argument));

            SPE * spe(0);
            {
                Lock l(*imp->mutex);

                imp->kernel_loaded->signal_and_wait(imp->mutex);
                spe = imp->spe;
            }
            std::cout << "SPEKernel " << imp->spe->id()<< ": loaded-signal received!" << std::endl;

            SPEEvent event(*imp->spe, SPE_EVENT_OUT_INTR_MBOX | SPE_EVENT_SPE_STOPPED);
            Instruction current_instruction;

            do
            {
                std::cout << "SPEKernel " << imp->spe->id() << ": Waiting for event!" << std::endl;
                spe_event_unit_t * e(event.wait(-1));
                std::cout << "SPEKernel " << imp->spe->id() << ": Event occured!" << std::endl;
                ASSERT(! (e->events & ~(SPE_EVENT_OUT_INTR_MBOX | SPE_EVENT_SPE_STOPPED)), "unexpected event happened!");

                if (e->events & SPE_EVENT_OUT_INTR_MBOX)
                {
                    std::cout << "SPEKernel " << imp->spe->id() << ": mail pending" << std::endl;
                    unsigned mail(0);

                    int retval(spe_out_intr_mbox_read(spe->context(), &mail, 1, SPE_MBOX_ALL_BLOCKING));
                    ASSERT(1 == retval, "weird return value in spe_*_mbox_read!");

                    do
                    {
                        switch (mail)
                        {
                            case km_result_dword:
                                std::cout << "SPEKernel " << imp->spe->id() << ": km_result_dword received!" << std::endl;
                                {
                                    Lock ll(*imp->mutex);

                                    spe_out_mbox_read(spe->context(),
                                            reinterpret_cast<unsigned int *>(imp->instructions[imp->spe_instruction_index].c.ea), 1);
                                }
                                continue;

                            case km_result_qword:
                                {
                                    Lock ll(*imp->mutex);

                                    spe_out_mbox_read(imp->spe->context(),
                                            reinterpret_cast<unsigned int *>(imp->instructions[imp->spe_instruction_index].c.ea), 2);
                                }
                                continue;

                            case km_instruction_finished:
                                std::cout << "SPEKernel " << imp->spe->id() << ": km_instruction_finished received!" << std::endl;
                                {
                                    Lock ll(*imp->mutex);

                                    imp->instructions[imp->spe_instruction_index].opcode = oc_noop;
                                    ++imp->finished_counter;
                                    ++imp->spe_instruction_index;
                                    imp->spe_instruction_index %= 8; /// \todo remove hardcoded numbers
                                    current_instruction = imp->instructions[imp->spe_instruction_index];
                                    imp->instruction_finished->broadcast();
                                }
                                continue;
                            case km_unknown_opcode:
                                std::cout << "SPEKernel " << imp->spe->id() << ": invalid opcode!" << std::endl;
                                throw std::string("Eek");
                                continue;
                        }

                        throw InternalError("Unexpected mail received in interrupt mailbox!");
                    } while (false);
                }
                else if (e->events & SPE_EVENT_SPE_STOPPED)
                {
                    std::cout << "SPEKernel: " << imp->spe->id() << " SPE stopped" << std::endl;
                    Lock ll(*imp->mutex);

                    imp->instruction_finished->broadcast();
                    break;
                }
                else
                {
                    std::cout << "SPEKernel " << imp->spe->id() << ": neither mail event nor stop event!" << std::endl;
                }
            } while (true);

            pthread_exit(0);
        }

        Implementation(const spe_program_handle_t & h, Environment * const e) :
            handle(h),
            environment(e),
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
            spe(NULL)
        {
            int retval(0);

            if (0 != posix_memalign(reinterpret_cast<void **>(&instructions), 128, 8 * sizeof(Instruction))) /// \todo remove hardcoded numbers
                throw std::bad_alloc(); /// \todo exception hierarchy

            if (0 != (retval = pthread_attr_init(attr)))
                throw ExternalError("libpthread", "pthread_attr_init failed, " + stringify(strerror(retval)));

            if (0 != (retval = pthread_attr_setdetachstate(attr, PTHREAD_CREATE_JOINABLE)))
                throw ExternalError("libpthread", "pthread_attr_setdetachstate failed," + stringify(strerror(retval)));

            if (0 != (retval = pthread_create(thread, 0, &kernel_thread, this)))
                throw ExternalError("libpthread", "pthread_create failed, " + stringify(strerror(retval)));

            Instruction noop = { oc_noop };
            std::fill(instructions, instructions + 8, noop);
        }

        ~Implementation()
        {
            int retval(0);

            if (0 != (retval = pthread_join(*thread, 0)))
                throw ExternalError("libpthread", "pthread_join failed, " + stringify(strerror(retval)));

            free(instructions);

            if (0 != (retval = pthread_attr_destroy(attr)))
                throw ExternalError("libpthread", "pthread_attr_destroy failed, " + stringify(strerror(retval)));

            delete attr;
            delete thread;
            delete instruction_finished;
            delete kernel_loaded;
            delete mutex;
        }

    };

    SPEKernel::SPEKernel(const spe_program_handle_t & handle, Environment * const environment) :
        _imp(new Implementation(handle, environment))
    {
    }

    SPEKernel::Iterator
    SPEKernel::begin() const
    {
        return Iterator(_imp->instructions);
    }

    SPEKernel::Iterator
    SPEKernel::end() const
    {
        return Iterator(_imp->instructions + 8); /// \todo remove hardcoded numbers
    }

    unsigned
    SPEKernel::enqueue(const Instruction & instruction)
    {
        CONTEXT("When enqueueing instruction:");
        unsigned result;
        Lock l(*_imp->mutex);

        while (_imp->spe_instruction_index == (_imp->next_free_index + 1) % 8) /// \todo remove hardcoded numbers
        {
            _imp->instruction_finished->wait(*_imp->mutex);
        }

        _imp->instructions[_imp->next_free_index] = instruction;
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
            throw ExternalError("libspe2", "spe_program_load failed, " + stringify(strerror(errno)));
        }

        _imp->spe = new SPE(spe);
        _imp->kernel_loaded->signal_and_wait(_imp->mutex);
        std::cout << "SPEKernel " << _imp->spe->id() << " : loading kernel into SPE" << std::endl;
        if (_imp->initial_instructions)
        {
            _imp->spe->signal();
        }
    }

    unsigned 
    SPEKernel::instruction_load() const
    {
        Lock l (*_imp->mutex);
        return (_imp->next_free_index - _imp->spe_instruction_index) % 8; /// \todo remove hardcoded numbers
    }

    void * SPEKernel::argument() const
    {
        Lock l(*_imp->mutex);

        return _imp->instructions;
    }

    void * const SPEKernel::environment() const
    {
        Lock l(*_imp->mutex);

        return _imp->environment;
    }
}
