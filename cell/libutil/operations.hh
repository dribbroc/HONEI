/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Util C++ library. LibUtil is free software;
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

#ifndef CELL_GUARD_UTIL_OPERATIONS_HH
#define CELL_GUARD_UTIL_OPERATIONS_HH 1

#include <cell/interface.hh>

namespace honei
{
    namespace cell
    {
        /**
         * \name Operation Data
         * \{
         *
         * Operation describes an operation that acts on a given set of entities.
         *
         * \ingroup grpframework
         */

        template <unsigned param_count_, typename DT_, ResultTransferMethod method_> struct Operation
        {
        };

        template <> struct Operation<1, float, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on a given block of data.
            void (*calculate)(vector float * elements, const unsigned size);
        };

        template <> struct Operation<1, float, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector float (*init)(void);

            /// Performs calculation on a given block of data.
            vector float (*calculate)(const vector float & accumulator, vector float * elements, const unsigned size);

            /// Finishs the calculations and returns a scalar result suitable for mail transfer.
            float (*finish)(const vector float & accumulator);
        };

        template <> struct Operation<2, float, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector float (*init)();

            /// Performs calculation on given blocks of data.
            vector float (*calculate)(const vector float & accumulator, vector float & carry,
                    vector float * a_elements, vector float * b_elements, const unsigned size, const unsigned offset);

            /// Finishs the calculations and returns a sclar result suitable for mail transfer.
            float (*finish)(const vector float & accumulator);
        };

        template <> struct Operation<2, float, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector float * a_elements, const vector float * b_elements, const unsigned size,
                    vector float & b_carry, const unsigned b_offset);
        };

        /// \}

        /**
         * \name Operation Skeletons
         * \{
         *
         * Function operation is a skeleton that executed a given operation's
         * methods.
         *
         * \ingroup grpframework
         */

        template <typename Operation_> typename Operation_::ResultType operation(const Operation_ & operation,
                const Instruction & instruction);

        template <> void operation<Operation<1, float, rtm_dma> >(const Operation<1, float, rtm_dma> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<1, float, rtm_mail> >(const Operation<1, float, rtm_mail> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<2, float, rtm_mail> >(const Operation<2, float, rtm_mail> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<2, float, rtm_dma> >(const Operation<2, float, rtm_dma> & operation,
                const Instruction & instruction);

        /// \}

        /// \name Utility functions often used by operations
        /// \{

        /**
         * Return a vector of biggest-possible floats.
         *
         * \ingroup grpframework
         */
        vector float biggest_float();

        /**
         * Return a vector of smalles-possible floats.
         *
         * \ingroup grpframework
         */
        vector float smallest_float();

        /**
         * Return a zero vector.
         *
         * \ingroup grpframework
         */
        vector float zero_float();

        /**
         * Return the sum of all of the vector's elements.
         *
         * \ingroup grpframework
         */
        float sum_float(const vector float & value);

        /**
         * Return the maximum of all of the vector's elements.
         *
         * \ingroup grpframework
         */
        float max_float(const vector float & value);

        /**
         * Return the minimum of all of the vector's elements.
         *
         * \ingroup grpframework
         */
        float min_float(const vector float & value);

        /// \}
    }
}

#endif
