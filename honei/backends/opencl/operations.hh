/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#ifndef OPENCL_GUARD_OPERATIONS_HH
#define OPENCL_GUARD_OPERATIONS_HH 1

#include <CL/cl.h>

namespace honei
{
    namespace opencl
    {
        void difference_float(void * r, void * x, void * y, unsigned long size, cl_device_type type);
        void difference_double(void * r, void * x, void * y, unsigned long size, cl_device_type type);

        float dot_product_float(void * x, void * y, unsigned long size, cl_device_type type);
        double dot_product_double(void * x, void * y, unsigned long size, cl_device_type type);

        void element_product_float(void * r, void * x, void * y, unsigned long size, cl_device_type type);
        void element_product_double(void * r, void * x, void * y, unsigned long size, cl_device_type type);

        void fill_float(void * x, float a, unsigned long size, cl_device_type type);
        void fill_double(void * x, double a, unsigned long size, cl_device_type type);

        float norm_l2_false_float(void * x, unsigned long size, cl_device_type type);
        double norm_l2_false_double(void * x, unsigned long size, cl_device_type type);

        void scale_float(void * x, float a, unsigned long size, cl_device_type type);
        void scale_double(void * x, double a, unsigned long size, cl_device_type type);

        void scaled_sum_float(void * r, void * x, void * y, float b, unsigned long size, cl_device_type type);
        void scaled_sum_double(void * r, void * x, void * y, double b, unsigned long size, cl_device_type type);
        void scaled_sum_float(void * r, void * x, void * y, unsigned long size, cl_device_type type);
        void scaled_sum_double(void * r, void * x, void * y, unsigned long size, cl_device_type type);

        void sum_float(void * r, void * x, void * y, unsigned long size, cl_device_type type);
        void sum_double(void * r, void * x, void * y, unsigned long size, cl_device_type type);

        void product_smell_dv_float(void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type);
        void product_smell_dv_double(void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type);

        void defect_smell_dv_float(void * rhs, void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type);
        void defect_smell_dv_double(void * rhs, void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type);
    }
}
#endif
