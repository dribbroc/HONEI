/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/attributes.hh>
#include <honei/cell/opcodes.hh>
#include <honei/cell/types.hh>

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
            km_instruction_finished = 1 << 0, ///< A completion message for an issued instruction.
            km_result_dword, ///< A result message carrying a double word of data.
            km_result_qword, ///< A result message carrying a quad word of data.
            km_unknown_opcode, ///< An error message about an unknown opcode.
            km_debug_get, ///< A debug message about a DMA GET transfer.
            km_debug_getl, ///< A debug message about a DMA GETL transfer.
            km_debug_put, ///< A debug message about a DMA PUT transfer.
            km_debug_putl, ///< A debug message about a DMA PUTL transfer.
            km_debug_enter, ///< A debug message about entering a function.
            km_debug_leave, ///< A debug message about leaving a function.
            km_debug_acquire_block, ///< A debug message about acquiring an LS memory block.
            km_debug_release_block, ///< A debug message about releasing an LS memory block.
            km_debug_dump, ///< A debug message to trigger a complete dump of the SPE's LS memory.
            km_debug_value, ///< A debug message for the output of a single value.
            km_profiler, ///< A profiler message carrying address and clock cycles for a given session.
            km_last = km_profiler
        };

        /**
         * KernelType identifies the operation type of a given SPE kernel.
         *
         * \ingroup grpspeinterface
         */
        enum KernelType
        {
            kt_stand_alone = 1 << 0, ///< A kernel that is not interconnected with other SPE kernels.
            kt_unknown, ///< An unknown (and thus invalid) kernel type.
            kt_last = kt_unknown
        };

        /**
         * Environment is used to inform a running SPE kernel about the begin and end of
         * its available data space in Local Store memory.
         *
         * \ingroup grpspeinterface
         */
        struct HONEI_ATTRIBUTE(packed) Environment
        {
            LocalStoreAddress begin; ///< The begin of the available data space in Local Store memory.

            LocalStoreAddress end; ///< The end of the available data space in Local Store memory.
        };

        /**
         * Operand is the data type of any operand specified for SPE instructions.
         *
         * \ingroup grpspeinterface
         */
        union HONEI_ATTRIBUTE(packed) Operand
        {
            EffectiveAddress ea; ///< An effective address pointing to a PPE-side memory location.

            long long s; ///< A signed 64 bit integer value.

            unsigned long long u; ///< An unsigned 64 bit integer value.

            double d; ///< A double precision floating point value.

            float f; ///< A single precision floating point value.

            float fa[2]; ///< An array of two single precision floating point values.

            /**
             * \name Contructors
             * \{
             *
             * Constructor.
             */

            Operand()
            {
            }

            Operand(const EffectiveAddress & value) :
                ea(value)
            {
            }

#ifdef __PPU__
            Operand(const unsigned long long & value) :
                u(value)
            {
            }
#endif

            Operand(const float & first, const float & second)
            {
                fa[0] = first;
                fa[1] = second;
            }

            /// \}
        };

        /**
         * Instruction is the data type used to convey an operation's OpCode and its
         * operands to the SPE-side kernel function.
         *
         * \ingroup grpspeinterface
         */
        struct HONEI_ATTRIBUTE(packed) Instruction
        {
            OpCode opcode; ///< An instruction's OpCode.

            unsigned size; ///< An instruction's default size.

            /// \name An instruction's Operands a to o
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
            Operand p;

            /// \}
        };

        /**
         * Capabilities is the data type used to inform about an SPE-side kernel's
         * supported opcodes as well as the type of the kernel.
         *
         * \ingroup grpspeinterface
         */
        struct HONEI_ATTRIBUTE(packed) Capabilities
        {
            KernelType type;

            unsigned opcode_count;

            const OpCode * opcodes;
        };

        /**
         * ResultTransferMethod enumerates all possible means to return data from the SPE
         * to the caller.
         *
         * \ingroup grpspeinterface
         */
        enum ResultTransferMethod
        {
            rtm_none, ///< The result cannot be transfered.
            rtm_mail, ///< The result will be transfered via mail.
            rtm_dma, ///< The result will be transfered via DMA.
        };

#if defined(__PPU__) || defined(DOXYGEN)

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
            const Capabilities capabilities; ///< A kernel's capabilities.

            const Environment environment; ///< A kernel's environment.

            const Handle handle; ///< A kernel's libspe2 handle.

            const std::string name; /// A kernel's name.
        };

#endif
    }
}

#endif
