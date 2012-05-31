
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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

#pragma once
#ifndef CUDA_GUARD_OPERATIONS_HH
#define CUDA_GUARD_OPERATIONS_HH 1

#include <cuda_runtime.h>

extern "C"
{
    /////////////////////////// LA /////////////////////////
    void cuda_defect_q1_float(void * rhs, void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_defect_q1_double(void * rhs, void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_defect_smell_dv_float(void * rhs, void * result, void * Aj, void * Ax, void * Arl, void * b,
            unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream = 0 );

    void cuda_defect_smell_dv_double(void * rhs, void * result, void * Aj, void * Ax, void * Arl, void * b,
            unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream = 0 );

    void cuda_defect_csr_dv_float(void * rhs, void * result, void * Aj, void * Ax, void * Ar, void * b,
            unsigned long rows, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream = 0);

    void cuda_defect_csr_dv_double(void * rhs, void * result, void * Aj, void * Ax, void * Ar, void * b,
            unsigned long rows, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream = 0);

    void cuda_difference_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize, cudaStream_t stream = 0);
    void cuda_difference_two_double(void * a, const void * b, unsigned long size, unsigned long blocksize, cudaStream_t stream = 0);
    void cuda_difference_three_float(void * r, const void * a, const void * b, unsigned long size, unsigned long blocksize);
    void cuda_difference_three_double(void * r, const void * a, const void * b, unsigned long size, unsigned long blocksize);

    float cuda_dot_product_two_float(const void * a, const void *b, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);
    double cuda_dot_product_two_double(const void * a, const void *b, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);

    void cuda_element_inverse_one_float(void * x, unsigned long size, unsigned long blocksize);
    void cuda_element_inverse_one_double(void * x, unsigned long size, unsigned long blocksize);

    void cuda_element_product_three_float(void * r, void * a, const void * b, unsigned long size, unsigned long blocksize);
    void cuda_element_product_three_double(void * r, void * a, const void * b, unsigned long size, unsigned long blocksize);

    float cuda_norm_l2_one_float(const void * a, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);

    double cuda_norm_l2_one_double(const void * a, unsigned long size, unsigned long blocksize,
            unsigned long gridsize);

    void cuda_product_bmdv_q1_float(void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_product_bmdv_q1_double(void * ll, void * ld, void * lu,
            void * dl, void * dd, void *du,
            void * ul, void * ud, void *uu, void * x, void * y,
            unsigned long size, unsigned long blocksize, unsigned long m);

    void cuda_product_smell_dv_float(void * x, void * y, void * Aj, void * Ax, void * Arl,
            unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream = 0);

    void cuda_product_smell_dv_double(void * x, void * y, void * Aj, void * Ax, void * Arl,
            unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
            unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream = 0);

    void cuda_product_csr_dv_float(void * x, void * y, void * Aj, void * Ax, void * Ar,
            unsigned long row_start, unsigned long row_end, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream = 0);

    void cuda_product_csr_dv_double(void * x, void * y, void * Aj, void * Ax, void * Ar,
            unsigned long row_start, unsigned long row_end, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream = 0);

    void cuda_prolongation_float(void * fine, unsigned long size_fine, void * coarse, unsigned long size_coarse,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_prolongation_double(void * fine, unsigned long size_fine, void * coarse, unsigned long size_coarse,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_restriction_float(void * coarse, unsigned long size_coarse, void * fine, unsigned long size_fine,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_restriction_double(void * coarse, unsigned long size_coarse, void * fine, unsigned long size_fine,
            unsigned long * macroBorderMask, unsigned long blocksize);

    void cuda_scaled_sum_three_float_s(void * x, void * y, void * z, float s, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_three_double_s(void * x, void * y, void * z, double s, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_three_float(void * x, void * y, void * z, unsigned long size, unsigned long blocksize);
    void cuda_scaled_sum_three_double(void * x, void * y, void * z, unsigned long size, unsigned long blocksize);

    void cuda_scale_one_float(void * x, const float a, unsigned long size, unsigned long blocksize);
    void cuda_scale_one_double(void * x, const double a, unsigned long size, unsigned long blocksize);

    void cuda_sum_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);
    void cuda_sum_two_double(void * a, const void * b, unsigned long size, unsigned long blocksize, cudaStream_t stream = 0);


    /////////////////////////// LBM /////////////////////////
    void cuda_eq_dist_grid_float(unsigned long start, unsigned long end, void * u, void * v, void * h,
            void * distribution_x, void * distribution_y,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            float g, float e,
            unsigned long blocksize);

    void cuda_eq_dist_grid_double(unsigned long start, unsigned long end, void * u, void * v, void * h,
            void * distribution_x, void * distribution_y,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            double g, double e,
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

    void cuda_collide_stream_grid_double(unsigned long start, unsigned long end,
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
            double tau, unsigned long size,
            unsigned long blocksize);

    void cuda_force_grid_float(
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * h, void * b,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float distribution_x_1, float distribution_y_1,
            float distribution_x_2, float distribution_y_2,
            float distribution_x_3, float distribution_y_3,
            float distribution_x_4, float distribution_y_4,
            float distribution_x_5, float distribution_y_5,
            float distribution_x_6, float distribution_y_6,
            float distribution_x_7, float distribution_y_7,
            float distribution_x_8, float distribution_y_8,
            float g, float d_x, float d_y, float d_t,
            unsigned long size,
            unsigned long blocksize);

    void cuda_force_grid_double(
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * h, void * b,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            double distribution_x_1, double distribution_y_1,
            double distribution_x_2, double distribution_y_2,
            double distribution_x_3, double distribution_y_3,
            double distribution_x_4, double distribution_y_4,
            double distribution_x_5, double distribution_y_5,
            double distribution_x_6, double distribution_y_6,
            double distribution_x_7, double distribution_y_7,
            double distribution_x_8, double distribution_y_8,
            double g, double d_x, double d_y, double d_t,
            unsigned long size,
            unsigned long blocksize);

    void cuda_force_grid_float_2(
            void * h, void * u, void * v,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float distribution_x_1, float distribution_y_1,
            float distribution_x_2, float distribution_y_2,
            float distribution_x_3, float distribution_y_3,
            float distribution_x_4, float distribution_y_4,
            float distribution_x_5, float distribution_y_5,
            float distribution_x_6, float distribution_y_6,
            float distribution_x_7, float distribution_y_7,
            float distribution_x_8, float distribution_y_8,
            float g, float d_x, float d_y, float d_t, float manning,
            unsigned long size,
            unsigned long blocksize);

    void cuda_force_grid_double_2(
            void * h, void * u, void * v,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            double distribution_x_1, double distribution_y_1,
            double distribution_x_2, double distribution_y_2,
            double distribution_x_3, double distribution_y_3,
            double distribution_x_4, double distribution_y_4,
            double distribution_x_5, double distribution_y_5,
            double distribution_x_6, double distribution_y_6,
            double distribution_x_7, double distribution_y_7,
            double distribution_x_8, double distribution_y_8,
            double g, double d_x, double d_y, double d_t, double manning,
            unsigned long size,
            unsigned long blocksize);

    void cuda_up_vel_dir_grid_float(unsigned long start, unsigned long end,
            void * types, void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            unsigned long blocksize);

    void cuda_up_vel_dir_grid_double(unsigned long start, unsigned long end,
            void * types, void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            unsigned long blocksize);

    void cuda_extraction_grid_dry_float(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, float epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_dry_double(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, double epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_wet_float(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, float epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_wet_double(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h, void * u, void * v,
            void * distribution_x, void * distribution_y, double epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_pollutant_dry_float(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h_flow, void * h_poll, float epsilon,
            unsigned long blocksize);

    void cuda_extraction_grid_pollutant_dry_double(
            unsigned long start, unsigned long end,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * h_flow, void * h_poll, double epsilon,
            unsigned long blocksize);

    ///////////////////////////  /////////////////////////
    void cuda_thread_synchronize();


    ///////////////////////////LBMFSI/////////////////////////
    void cuda_collide_stream_fsi_float(unsigned long start, unsigned long end,
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            void * f_mea_1, void * f_mea_2,
            void * f_mea_3, void * f_mea_4, void * f_mea_5,
            void * f_mea_6, void * f_mea_7, void * f_mea_8,
            void * line_flags, void * dist_x, void * dist_y,
            float d_xu, float d_yv,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_1_float(
            void * cuda_dir_5,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_2_float(
            void * cuda_dir_6,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_3_float(
            void * cuda_dir_7,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_4_float(
            void * cuda_dir_8,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_5_float(
            void * cuda_dir_1,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_6_float(
            void * cuda_dir_2,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_7_float(
            void * cuda_dir_3,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);

    void cuda_boundary_init_fsi_dir_8_float(
            void * cuda_dir_4,
            void * boundary_flags ,
            void * solid_flags ,
            void * solid_old_flags ,
            void * h ,
            void * u ,
            void * v ,
            void * f_0 ,
            void * f_1 ,
            void * f_2 ,
            void * f_3 ,
            void * f_4 ,
            void * f_5 ,
            void * f_6 ,
            void * f_7 ,
            void * f_8 ,
            void * f_eq_0 ,
            void * f_eq_1 ,
            void * f_eq_2 ,
            void * f_eq_3 ,
            void * f_eq_4 ,
            void * f_eq_5 ,
            void * f_eq_6 ,
            void * f_eq_7 ,
            void * f_eq_8 ,
            void * f_temp_0 ,
            void * f_temp_1 ,
            void * f_temp_2 ,
            void * f_temp_3 ,
            void * f_temp_4 ,
            void * f_temp_5 ,
            void * f_temp_6 ,
            void * f_temp_7 ,
            void * f_temp_8 ,
            float current_u,
            float current_v,
            unsigned long size,
            unsigned long blocksize);
}
#endif
