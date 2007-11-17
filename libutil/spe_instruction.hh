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

namespace honei
{
    class SPEKernel;

    class SPEInstruction
    {
        private:
            /// Our kernel.
            mutable SPEKernel * _kernel;

            /// Our index.
            mutable unsigned _index;

            /// Our instruction structure.
            Instruction _instruction;

            static const Operand empty;

            /// Enqueue us with a given kernel.
            void enqueue_with(SPEKernel * kernel) const;

        public:
            friend class SPEManager::Implementation;

            /// Constructor.
            SPEInstruction() { }

            /// Constructor.
            SPEInstruction(const OpCode opcode, const unsigned size, const Operand a = empty,
                    const Operand b = empty, const Operand c = empty, const Operand d = empty,
                    const Operand e = empty, const Operand f = empty, const Operand g = empty,
                    const Operand h = empty, const Operand i = empty, const Operand j = empty,
                    const Operand k = empty);

            /// Destructor.
            ~SPEInstruction();

            /// Wait until we have been executed.
            void wait();

            /// Returns our instruction
            Instruction instruction()
            {
                return _instruction;
            }
    };
}

#endif
