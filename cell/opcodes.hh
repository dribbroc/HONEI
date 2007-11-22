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

#ifndef LIBLA_CELL_GUARD_OPCODES_HH
#define LIBLA_CELL_GUARD_OPCODES_HH 1

namespace honei
{
    namespace cell
    {
        enum OpCode
        {
            oc_noop = 1 << 0,
            oc_halt,

            oc_dense_dense_float_sum = 1 << 4,
            oc_dense_dense_float_difference,
            oc_dense_dense_float_dot_product,
            oc_dense_dense_float_matrix_product,
            oc_dense_dense_float_matrix_vector_product,
            oc_banded_dense_float_matrix_vector_product,
            oc_dense_dense_float_element_product,
            oc_dense_sparse_float_sum,
            oc_dense_float_reduction_sum,
            oc_dense_float_reduction_max,
            oc_dense_float_reduction_min,
            oc_dense_float_scale,

            oc_test_instruction_finished = 1 << 30,
            oc_test_result_dword,
            oc_test_result_qword,

            oc_last = oc_test_result_qword
        };
    }
}

#endif
