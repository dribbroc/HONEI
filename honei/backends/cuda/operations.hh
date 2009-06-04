
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CUDA_GUARD_OPERATIONS_HH
#define CUDA_GUARD_OPERATIONS_HH 1

extern "C"
{
    /////////////////////////// LA /////////////////////////
    void cuda_defect_q1_float(void * rhs, void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_defect_smell_dv_float(void * rhs, void * result, void * Aj, void * Ax, void * b,
            unsigned long rows, unsigned long columns, unsigned long num_cols_per_row, unsigned long stride,
            unsigned long blocksize);

    void cuda_defect_smell_dv_double(void * rhs, void * result, void * Aj, void * Ax, void * b,
            unsigned long rows, unsigned long columns, unsigned long num_cols_per_row, unsigned long stride,
            unsigned long blocksize);

    void cuda_difference_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);

    float cuda_dot_product_two_float(const void * a, const void *b, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);

    void cuda_element_inverse_one_float(void * x, unsigned long size, unsigned long blocksize);

    void cuda_element_product_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);

    float cuda_norm_l2_one_float(const void * a, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);

    void cuda_product_bmdv_q1_float(void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_product_smell_dv_float(void * x, void * y, void * Aj, void * Ax,
            unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize);

    void cuda_product_smell_dv_double(void * x, void * y, void * Aj, void * Ax,
            unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize);

    void cuda_prolongation_float(void * fine, unsigned long size_fine, void * coarse, unsigned long size_coarse,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_restriction_float(void * coarse, unsigned long size_coarse, void * fine, unsigned long size_fine,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_scaled_sum_two_float(void * x, const void * y, float b, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_two_double(void * x, const void * y, double b, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_three_float(void * x, void * y, void * z, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_three_double(void * x, void * y, void * z, unsigned long size, unsigned long blocksize);

    void cuda_scale_one_float(void * x, const float a, unsigned long size, unsigned long blocksize);

    void cuda_sum_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);


    /////////////////////////// LBM /////////////////////////
    void cuda_eq_dist_grid_float(unsigned long start, unsigned long end, void * u, void * v, void * h,
            void * distribution_x, void * distribution_y,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            float g, float e,
            unsigned long blocksize);

    void cuda_collide_stream_grid_float(unsigned long start, unsigned long end,
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * f_temp_0, void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float tau, unsigned long size,
            unsigned long blocksize);

    void cuda_force_grid_float(
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * h, void * b,
            void * distribution_x, void * distribution_y,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float g, float d_x, float d_y, float d_t,
            unsigned long size,
            unsigned long blocksize);

    void cuda_up_vel_dir_grid_float(void * types,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            unsigned long size,
            unsigned long blocksize);

    void cuda_extraction_grid_dry_float(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, float epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_wet_float(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, float epsilon,
            unsigned long blocksize);

    ///////////////////////////  /////////////////////////
    void cuda_thread_synchronize();
}
#endif
