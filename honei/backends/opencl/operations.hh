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

#ifndef OPENCL_GUARD_OPERATIONS_HH
#define OPENCL_GUARD_OPERATIONS_HH 1

#include <CL/cl.h>

namespace honei
{
    namespace opencl
    {
        void scaled_sum_float(void * x, void * y, float b, unsigned long size, cl_device_type type);
        void scaled_sum_double(void * x, void * y, double b, unsigned long size, cl_device_type type);

        void product_smell_dv_float(void * x, void * y, void * Aj, void * Ax,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, cl_device_type type);
        void product_smell_dv_double(void * x, void * y, void * Aj, void * Ax,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, cl_device_type type);

        void copy_float(void * x, void * y, unsigned long size, cl_device_type type);
        void copy_double(void * x, void * y, unsigned long size, cl_device_type type);
    }
}
#endif
