/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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

#ifndef LIBUTIL_GUARD_SPE_INSTRUCTION_HH
#define LIBUTIL_GUARD_SPE_INSTRUCTION_HH 1

#include <honei/cell/cell.hh>
#include <honei/libutil/spe_manager.hh>

#include <list>
#include <vector>

namespace honei
{
    class SPEKernel;
    class SPEInstructionQueue;
    class SPEInstructionStream;

    class SPEInstruction
    {
        public:
            typedef honei::cell::Instruction Instruction;
            typedef honei::cell::Operand Operand;
            typedef honei::cell::OpCode OpCode;

        protected:
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

            static const Operand empty;

            /// Enqueue us with a given kernel.
            void _enqueue_with(SPEKernel * kernel) const;

        public:
            friend class SPEManager::Implementation;
            friend class SPEInstructionQueue;
            friend class SPEInstructionStream;

            /// Constructor.
            SPEInstruction(const OpCode opcode, const unsigned size, const Operand & a = empty,
                    const Operand & b = empty, const Operand & c = empty, const Operand & d = empty,
                    const Operand & e = empty, const Operand & f = empty, const Operand & g = empty,
                    const Operand & h = empty, const Operand & i = empty, const Operand & j = empty,
                    const Operand & k = empty, const Operand & l = empty, const Operand & m = empty,
                    const Operand & n = empty, const Operand & o = empty);

            /// Destructor.
            ~SPEInstruction();

            /// Wait until we have been executed.
            void wait() const;

            /// Returns true if the instruction has been finished
            bool finished() const;

            /// Returns our instruction
            Instruction & instruction() const;
    };

    template <unsigned number_of_entities_, typename DataType_, cell::ResultTransferMethod method_> class SPEFrameworkInstruction;

    template <typename DataType_> class SPEFrameworkInstruction<1, DataType_, cell::rtm_dma> :
        public SPEInstruction
    {
        private:
            typedef cell::Instruction Instruction;

            typedef cell::OpCode OpCode;

            /// Index of the first element that is to be transfered.
            unsigned _begin_transfers;

            /// Index of the element behind the last element that is to be transfered.
            unsigned _end_transfers;

            /// Will the SPE be used?
            bool _use_spe;

        public:
            /**
             * Constructor.
             *
             * \param opcode The instruction's opcode.
             * \param elements The pointer to the elements of the container.
             * \param size The overall size of the container.
             */
            SPEFrameworkInstruction(const OpCode opcode, DataType_ * elements, const unsigned size);

            unsigned transfer_begin() const
            {
                return _begin_transfers;
            }

            unsigned transfer_end() const
            {
                return _end_transfers;
            }

            bool use_spe() const
            {
                return _use_spe;
            }
    };

    template <typename DataType_> class SPEFrameworkInstruction<1, DataType_, cell::rtm_mail> :
        public SPEInstruction
    {
        private:
            typedef cell::Instruction Instruction;

            typedef cell::OpCode OpCode;

            /// Index of the first element that is to be transfered.
            unsigned _begin_transfers;

            /// Index of the element behind the last element that is to be transfered.
            unsigned _end_transfers;

            /// Will the SPE be used?
            bool _use_spe;

        public:
            /**
             * Constructor.
             *
             * \param opcode The instruction's opcode.
             * \param result The pointer to the pre-allocated result-memory.
             * \param elements The pointer to the elements of the container.
             * \param size The overall size of the container.
             */
            SPEFrameworkInstruction(const OpCode opcode, DataType_ * result, DataType_ * elements, const unsigned size);

            unsigned transfer_begin() const
            {
                return _begin_transfers;
            }

            unsigned transfer_end() const
            {
                return _end_transfers;
            }

            bool use_spe() const
            {
                return _use_spe;
            }
    };

    extern template class SPEFrameworkInstruction<1, float, cell::rtm_dma>;

    extern template class SPEFrameworkInstruction<1, float, cell::rtm_mail>;

    template <typename DataType_> class SPEFrameworkInstruction<2, DataType_, cell::rtm_dma> :
        public SPEInstruction
    {
        private:
            typedef cell::Instruction Instruction;

            typedef cell::OpCode OpCode;

            /// Index of the first element that is to be transfered.
            unsigned _begin_transfers;

            /// Index of the element behind the last element that is to be transfered.
            unsigned _end_transfers;

            /// Will the SPE be used?
            bool _use_spe;

        public:
            /**
             * Constructor.
             *
             * \param opcode The instruction's opcode.
             * \param result The pointer to the pre-allocated result-memory.
             * \param elements The pointer to the elements of the container.
             * \param size The overall size of the container.
             */
            SPEFrameworkInstruction(const OpCode opcode, DataType_ * result, DataType_ * elements, const unsigned size);

            unsigned transfer_begin() const
            {
                return _begin_transfers;
            }

            unsigned transfer_end() const
            {
                return _end_transfers;
            }

            bool use_spe() const
            {
                return _use_spe;
            }
    };

    extern template class SPEFrameworkInstruction<2, float, cell::rtm_dma>;

    class SPEInstructionQueue
    {
        private:
            struct Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// Enqueue us with a given kernel.
            void _enqueue_with(SPEKernel * kernel);

        public:
            friend class SPEManager::Implementation;

            /// Constructor.
            SPEInstructionQueue();

            /// Destructor.
            ~SPEInstructionQueue();

            /// Insert an instruction at the end of the queue.
            void push_back(const SPEInstruction & instruction);

            /// Wait until all our instructions have been executed.
            void wait() const;

            /// Returns an iterator pointing to the front of our queue
            std::list<SPEInstruction>::iterator begin();

            /// Returns an iterator pointing behind the last element of our queue
            std::list<SPEInstruction>::iterator end();

            /// Returns our instruction count.
            unsigned long size();

            /// Returns true if all instructions have been finished
            bool finished() const;
    };

    class SPEInstructionStream
    {
        private:
            struct Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// Enqueue us with a given kernel.
            void _enqueue_with(SPEKernel * kernel);

        public:
            friend class SPEManager::Implementation;

            /// Constructor.
            SPEInstructionStream(const SPEInstruction & instruction);

            /// Destructor.
            ~SPEInstructionStream();

            /// Insert an instruction at the end of the queue.
            void input(const SPEInstruction & instruction);

            /// Wait until all our instructions have been executed.
            void wait() const;

            /// Returns an iterator pointing to the front of our stream
            std::vector<SPEInstruction>::iterator begin();

            /// Returns an iterator pointing behind the last element of our stream
            std::vector<SPEInstruction>::iterator end();

            /// Returns our instruction count.
            unsigned long size();

            /// Returns true if all instructions have been finished
            bool finished() const;
    };
}
#endif
