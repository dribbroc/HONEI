/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CELL_GUARD_INTERFACE_HH
#define CELL_GUARD_INTERFACE_HH 1

#include <cell/opcodes.hh>
#include <cell/types.hh>

#ifdef __PPU__
#  include <libspe2-types.h>
#  include <string>
#endif

namespace honei
{
    namespace cell
    {
        /**
         * KernelMessages is the set of unique ids to identify a kernel message
         * from SPE to PPE.
         *
         * \ingroup grpspeinterface
         */
        enum KernelMessages
        {
            km_instruction_finished = 1 << 0, /// < A completion message for an issued instruction.
            km_result_dword, /// < A result message carrying a double word of data.
            km_result_qword, /// < A result message carrying a quad word of data.
            km_unknown_opcode, /// < An error message about an unknown opcode.
            km_debug_get, /// < A debug message about a DMA GET transfer.
            km_debug_put, /// < A debug message about a DMA PUT transfer.
            km_debug_enter, /// < A debug message about entering a function.
            km_debug_leave, /// < A debug message about leaving a function.
            km_debug_acquire_block, /// < A debug message about acquiring an LS memory block.
            km_debug_release_block, /// < A debug message about releasing an LS memory block.
            km_debug_dump, /// < A debug message to trigger a complete dump of the SPE's LS memory.
            km_last = km_debug_dump
        };

        /**
         * KernelType identifies the operation type of a given SPE kernel.
         *
         * \ingroup grpspeinterface
         */
        enum KernelType
        {
            kt_stand_alone = 1 << 0, /// < A kernel that is not interconnected with other SPE kernels.
            kt_unknown, /// < An unknown (and thus invalid) kernel type.
            kt_last = kt_unknown
        };

        /**
         * Environment is used to a running SPE kernel about the begin and end of
         * its available data space in Local Store memory.
         *
         * \ingroup grpspeinterface
         */
        struct __attribute__((packed)) Environment
        {
            LocalStoreAddress begin; /// < The begin of the available data space in Local Store memory.

            LocalStoreAddress end; /// < The end of the available data space in Local Store memory.
        };

        /**
         * Operand is the data type of any operand specified for SPE instructions.
         *
         * \ingroup grpspeinterface
         */
        union Operand
        {
            EffectiveAddress ea; /// < An effective address pointing to a PPE-side memory location.

            long long s; /// < A signed 64 bit integer value.

            unsigned long long u; /// < An unsigned 64 bit integer value.

            double d; /// < A double precision floating point value.

            float f; /// < A double precision floating point value.
        };

        /**
         * Instruction is the data type used to convey an operation's OpCode and its
         * operands to the SPE-side kernel function.
         *
         * \ingroup grpspeinterface
         */
        struct __attribute__((packed)) Instruction
        {
            OpCode opcode; /// < An instruction's OpCode.

            unsigned size; /// < An instruction's default size.

            /// \name An instructions's Operands a to o
            /// \{

            Operand a;
            Operand b;
            Operand c;
            Operand d;
            Operand e;
            Operand f;
            Operand g;
            Operand h;
            Operand i;
            Operand j;
            Operand k;
            Operand l;
            Operand m;
            Operand n;
            Operand o;

            /// \}
        };

        /**
         * Capabilities is the data type used to inform about an SPE-side kernel's
         * supported opcodes as well as the type of the kernel.
         *
         * \ingroup grpspeinterface
         */
        struct __attribute__((packed)) Capabilities
        {
            KernelType type;

            unsigned opcode_count;

            const OpCode * opcodes;
        };

#ifdef __PPU__

        /**
         * Handle is the data type used to load a kernel program into an SPE.
         *
         * \ingroup grpspeinterface
         */
        typedef const ::spe_program_handle_t & Handle;

        /**
         * KernelInfo is the data type used to gather all information to create a new
         * kernel.
         *
         * \ingroup grpspeinterface
         */
        struct KernelInfo
        {
            const Capabilities capabilities; /// < A kernel's capabilities.

            const Environment environment; /// < A kernel's environment.

            const Handle handle; /// < A kernel's libspe2 handle.

            const std::string name; /// A kernel's name.
        };

#endif
    }
}

#endif
