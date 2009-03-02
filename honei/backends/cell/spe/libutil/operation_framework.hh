/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef CELL_GUARD_UTIL_OPERATION_FRAMEWORK_HH
#define CELL_GUARD_UTIL_OPERATION_FRAMEWORK_HH 1

#include <honei/backends/cell/interface.hh>

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
            void (*calculate)(vector float * elements, const unsigned size, const float optional_scalar);
        };

        template <> struct Operation<1, double, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on a given block of data.
            void (*calculate)(vector double * elements, const unsigned size, const double optional_scalar);
        };

        template <> struct Operation<1, float, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector float (*init)(void);

            /// Performs calculation on a given block of data.
            vector float (*calculate)(const vector float & accumulator, vector float * elements, const unsigned size,
                    const float optional_scalar);

            /// Finishs the calculations and returns a scalar result suitable for mail transfer.
            float (*finish)(const vector float & accumulator);
        };

        template <> struct Operation<1, double, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector double (*init)(void);

            /// Performs calculation on a given block of data.
            vector double (*calculate)(const vector double & accumulator, vector double * elements, const unsigned size,
                    const double optional_scalar);

            /// Finishs the calculations and returns a scalar result suitable for mail transfer.
            double (*finish)(const vector double & accumulator);
        };

        template <> struct Operation<2, float, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector float * a_elements, const vector float * b_elements, const unsigned size,
                    vector float & b_carry, const unsigned b_offset, const float optional_scalar);
        };

        template <> struct Operation<2, double, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector double * a_elements, const vector double * b_elements, const unsigned size,
                    vector double & b_carry, const unsigned b_offset, const double optional_scalar);
        };

        template <> struct Operation<2, float, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector float (*init)();

            /// Performs calculation on given blocks of data.
            vector float (*calculate)(const vector float & accumulator, vector float & carry, const vector float * a_elements,
                    const vector float * b_elements, unsigned size, unsigned offset, float scalar1, float scalar2);

            /// Finishs the calculations and returns a sclar result suitable for mail transfer.
            float (*finish)(const vector float & accumulator);
        };

        template <> struct Operation<2, double, rtm_mail>
        {
            /// Our result type.
            typedef unsigned ResultType;

            /// Initialises the accumulator.
            vector double (*init)();

            /// Performs calculation on given blocks of data.
            vector double (*calculate)(const vector double & accumulator, vector double & carry, const vector double * a_elements,
                    const vector double * b_elements, unsigned size, unsigned offset, double scalar1, double scalar2);

            /// Finishs the calculations and returns a sclar result suitable for mail transfer.
            double (*finish)(const vector double & accumulator);
        };

        template <> struct Operation<3, float, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector float * a_elements, const vector float * b_elements, const vector float * c_elements,
                    const unsigned size, vector float & b_carry, const unsigned b_offset, vector float & c_carry,
                    const unsigned c_offset, const float optional_scalar);
        };

        template <> struct Operation<3, double, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector double * a_elements, const vector double * b_elements, const vector double * c_elements,
                    const unsigned size, vector double & b_carry, const unsigned b_offset, vector double & c_carry,
                    const unsigned c_offset, const double optional_scalar);
        };

        template <> struct Operation<4, float, rtm_dma>
        {
            /// Our result type.
            typedef void ResultType;

            /// Performs initialization.
            void (*init)(float scalar, vector float carry);

            /// Performs calculation on given blocks of data.
            void (*calculate)(vector float * result_elements, const vector float * a_elements, const vector float * b_elements,
                    const vector float * c_elements, const unsigned size, float scalar);
        };

        /// \}

        /**
         * \name Operation Skeletons
         * \{
         *
         * Function operation is a skeleton that executes a given operation's
         * methods.
         *
         * \ingroup grpframework
         */

        template <typename Operation_> typename Operation_::ResultType operation(const Operation_ & operation,
                const Instruction & instruction);

        template <> void operation<Operation<1, float, rtm_dma> >(const Operation<1, float, rtm_dma> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<1, double, rtm_dma> >(const Operation<1, double, rtm_dma> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<1, float, rtm_mail> >(const Operation<1, float, rtm_mail> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<1, double, rtm_mail> >(const Operation<1, double, rtm_mail> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<2, float, rtm_dma> >(const Operation<2, float, rtm_dma> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<2, double, rtm_dma> >(const Operation<2, double, rtm_dma> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<2, float, rtm_mail> >(const Operation<2, float, rtm_mail> & operation,
                const Instruction & instruction);

        template <> unsigned operation<Operation<2, double, rtm_mail> >(const Operation<2, double, rtm_mail> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<3, float, rtm_dma> >(const Operation<3, float, rtm_dma> & operation,
                const Instruction & instruction);

        template <> void operation<Operation<4, float, rtm_dma> >(const Operation<4, float, rtm_dma> & operation,
                const Instruction & instruction);
        /// \}

        /// \name Utility functions often used by operations
        /// \{

        /**
         * \{
         * Return a vector of biggest-possible values
         *
         * \ingroup grpframework
         */

        vector float biggest_float();

        vector double biggest_double();

        /// \}

        /**
         * \{
         * Return a vector of smallesit-possible values
         *
         * \ingroup grpframework
         */

        vector float smallest_float();

        vector double smallest_double();

        /// \}

        /**
         * \{
         * Return a zero vector.
         *
         * \ingroup grpframework
         */

        vector float zero_float();

        vector double zero_double();

        /// \}

        /**
         * \{
         * Return the sum of all of the vector's elements.
         *
         * \ingroup grpframework
         */

        float sum_float(const vector float & value);

        double sum_double(const vector double & value);

        /// \}

        /**
         * \{
         * Return the maximum of all of the vector's elements.
         *
         * \ingroup grpframework
         */

        float max_float(const vector float & value);

        double max_double(const vector double & value);

        /// \}

        /**
         * \{
         * Return the minimum of all of the vector's elements.
         *
         * \ingroup grpframework
         */

        float min_float(const vector float & value);

        double min_double(const vector double & value);

        /// \}

        /// \}
    }
}

#endif
