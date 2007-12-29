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

#ifndef LIBUTIL_GUARD_SPE_INSTRUCTION_HH
#define LIBUTIL_GUARD_SPE_INSTRUCTION_HH 1

#include <cell/cell.hh>
#include <libutil/spe_manager.hh>

#include <list>

namespace honei
{
    class SPEKernel;
    class SPEInstructionQueue;

    class SPEInstruction
    {
        public:
            typedef honei::cell::Instruction Instruction;
            typedef honei::cell::Operand Operand;
            typedef honei::cell::OpCode OpCode;

        private:
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

            static const Operand empty;

            /// Enqueue us with a given kernel.
            void _enqueue_with(SPEKernel * kernel) const;

        public:
            friend class SPEManager::Implementation;
            friend class SPEInstructionQueue;

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

            /// Returns our instruction
            Instruction instruction() const;
    };

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

            /// Insert an instruction at the end of the queue.
            const SPEInstruction & push_back(const SPEInstruction & instruction);

            /// Wait until all our instructions have been executed.
            void wait() const;

            /// Returns an iterator pointing to the front of our queue
            std::list<SPEInstruction>::iterator begin();

            /// Returns an iterator pointing to the front of our queue
            std::list<SPEInstruction>::iterator end();

            /// Returns our instruction count.
            unsigned long size();
    };
}
#endif
