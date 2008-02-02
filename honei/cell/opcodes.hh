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

            // LibLA opcodes
            oc_sum_dense_dense_float = 1 << 4,
            oc_sum_dense_dense_double,
            oc_sum_dense_scalar_float,
            oc_sum_dense_scalar_double,
            oc_difference_dense_dense_float,
            oc_difference_dense_dense_double,
            oc_dot_product_dense_dense_float,
            oc_dot_product_dense_dense_double,
            oc_product_dense_matrix_dense_matrix_float,
            oc_product_dense_matrix_dense_vector_float,
            oc_product_banded_matrix_dense_vector_float,
            oc_product_banded_matrix_dense_vector_double,
            oc_element_product_dense_dense_float,
            oc_element_product_dense_dense_double,
            oc_sum_dense_sparse_float,
            oc_reduction_sum_dense_float,
            oc_reduction_sum_dense_double,
            oc_reduction_max_dense_float,
            oc_reduction_min_dense_float,
            oc_scale_dense_float,
            oc_norm_max_dense_float,
            oc_norm_l_one_dense_float,
            oc_norm_l_two_dense_float,
            oc_element_inverse_float,
            oc_scaled_sum_dense_dense_float,
            oc_scale_sparse_float,
            oc_fill_float,
            oc_fill_double,

            // LibGraph opcodes
            oc_node_distance_float_wkk,
            oc_node_distance_float_wfr,

            // LibMath opcodes
            oc_sqrt_dense_float,

            // LibSWE opcodes
            oc_flow_processing_x_float,
            oc_flow_processing_y_float,

            // Unittest opcodes
            oc_test_instruction_finished = 1 << 30,
            oc_test_result_dword,
            oc_test_result_qword,

            oc_last = oc_test_result_qword
        };
    }
}

#endif
