/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef SSE_GUARD_OPERATIONS_HH
#define SSE_GUARD_OPERATIONS_HH 1

namespace honei
{
    namespace sse
    {
        ///////////// LA
        void difference(float * a, const float * b, unsigned long size);
        void difference(double * a, const double * b, unsigned long size);

        float dot_product(const float * a, float * b, unsigned long size);
        double dot_product(double * a, double * b, unsigned long size);

        void element_inverse(float * x, unsigned long size);
        void element_inverse(double * x, unsigned long size);

        void element_product(float * a, const float * b, unsigned long size);
        void element_product(double * a, const double * b, unsigned long size);
        void element_product(float * a, const float * b, const float * c, unsigned long size);
        void element_product(double * a, const double * b, const double * c, unsigned long size);

        float norm_l2(const float * a, unsigned long size);
        double norm_l2(double * a, unsigned long size);

        void product_dm(float * x, float * y, float b, unsigned long size);
        void product_dm(double * x, double * y, double b, unsigned long size);
        void product_dm_nx2(float * result, const float * a, const float * b, unsigned long size);
        void product_dm_nx2(double * result, const double * a, const double * b, unsigned long size);
        void product_bmdv(float * x, const float * y, const float * z, unsigned long size);
        void product_bmdv(double * x, const double * y, const double * z, unsigned long size);
        void product_bmdv_q1(float * ll, float * ld, float * lu,
                float * dl, float * dd, float * du,
                float * ul, float * ud, float * uu,
                float * b, float * result,
                unsigned long, unsigned long m);
        void product_bmdv_q1(double * ll, double * ld, double * lu,
                double * dl, double * dd, double * du,
                double * ul, double * ud, double * uu,
                double * b, double * result,
                unsigned long, unsigned long m);
        void product_smell_dv(float * result, unsigned long * Aj, float * Ax, float * b,
            unsigned long stride, unsigned long rows, unsigned long num_cols_per_row);
        void product_smell_dv(float * result, unsigned long * Aj, float * Ax, unsigned long * Arl, float * b,
            unsigned long stride, unsigned long rows, unsigned long num_cols_per_row,
            unsigned long row_start, unsigned long row_end);
        void product_smell_dv(double * result, unsigned long * Aj, double * Ax, double * b,
            unsigned long stride, unsigned long rows, unsigned long num_cols_per_row);
        void product_smell_dv(double * result, unsigned long * Aj, double * Ax, unsigned long * Arl, double * b,
            unsigned long stride, unsigned long rows, unsigned long num_cols_per_row,
            unsigned long row_start, unsigned long row_end);

        float reduction_sum(const float * a, unsigned long size);
        double reduction_sum(double * a, unsigned long size);

        void scaled_sum(float * x, const float * y, float b, unsigned long size);
        void scaled_sum(double * x, const double * y, double b, unsigned long size);
        void scaled_sum(float * x, const float * y, const float * z, unsigned long size);
        void scaled_sum(double * x, const double * y, const double * z, unsigned long size);
        void scaled_sum(float * x, const float * y, const float * z, float b, unsigned long size);
        void scaled_sum(double * x, const double * y, const double * z, double b, unsigned long size);

        void scale(const float a, float * x, unsigned long size);
        void scale(const double a, double * x, unsigned long size);

        void sum(float * a, const float * b, unsigned long size);
        void sum(double * a, const double * b, unsigned long size);
        void sum(const float a, float * x, unsigned long size);
        void sum(const double a, double * x, unsigned long size);

        ///////////// LBM
        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * f_eq_0);

        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * distribution_x, float * distribution_y,
                float * f_eq,
                unsigned long dir);

        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * distribution_x, float * distribution_y,
                float * f_eq,
                unsigned long dir);

        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * f_eq_0);

        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * distribution_x, double * distribution_y,
                double * f_eq,
                unsigned long dir);

        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * distribution_x, double * distribution_y,
                double * f_eq,
                unsigned long dir);

        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, float tau,
                float * f_temp_0, float * f_0, float * f_eq_0);

        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, double tau,
                double * f_temp_0, double * f_0, double * f_eq_0);

        void collide_stream_grid_dir_n(unsigned long end, float tau,
                unsigned long * dir, unsigned long * dir_index,
                float * f_temp, float * f, float * f_eq);

        void collide_stream_grid_dir_n(unsigned long end, double tau,
                unsigned long * dir, unsigned long * dir_index,
                double * f_temp, double * f, double * f_eq);

        void extraction_grid_dry(unsigned long begin, unsigned long end,
                float * distribution_x, float * distribution_y,
                float * h, float * u, float * v,
                float * f_0, float * f_1, float * f_2,
                float * f_3, float * f_4, float * f_5,
                float * f_6, float * f_7, float * f_8, float epsilon);

        void extraction_grid_dry(unsigned long begin, unsigned long end,
                double * distribution_x, double * distribution_y,
                double * h, double * u, double * v,
                double * f_0, double * f_1, double * f_2,
                double * f_3, double * f_4, double * f_5,
                double * f_6, double * f_7, double * f_8, double epsilon);

        void extraction_grid_wet(unsigned long begin, unsigned long end,
                float * distribution_x, float * distribution_y,
                float * h, float * u, float * v,
                float * f_0, float * f_1, float * f_2,
                float * f_3, float * f_4, float * f_5,
                float * f_6, float * f_7, float * f_8, float epsilon);

        void extraction_grid_wet(unsigned long begin, unsigned long end,
                double * distribution_x, double * distribution_y,
                double * h, double * u, double * v,
                double * f_0, double * f_1, double * f_2,
                double * f_3, double * f_4, double * f_5,
                double * f_6, double * f_7, double * f_8, double epsilon);

        void up_vel_dir_grid(unsigned long begin, unsigned long end, unsigned long * limits, unsigned long * types,
                float * f_temp_1, float * f_temp_2,  float * f_temp_3,
                float * f_temp_4, float * f_temp_5,  float * f_temp_6,
                float * f_temp_7, float * f_temp_8);

        void up_vel_dir_grid(unsigned long begin, unsigned long end, unsigned long * limits, unsigned long * types,
                double * f_temp_1, double * f_temp_2,  double * f_temp_3,
                double * f_temp_4, double * f_temp_5,  double * f_temp_6,
                double * f_temp_7, double * f_temp_8);
    }
}

#endif
