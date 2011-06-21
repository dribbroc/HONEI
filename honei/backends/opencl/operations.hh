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
#include <string>
#include <honei/backends/opencl/opencl_backend.hh>

namespace honei
{
    namespace opencl
    {
        void difference(void * r, void * x, void * y, unsigned long size, cl_device_type type, std::string function);

        template <typename DT_>
        DT_ dot_product(void * x, void * y, unsigned long size, cl_device_type type, std::string function);

        void element_product(void * r, void * x, void * y, unsigned long size, cl_device_type type, std::string function);

        void fill_float(void * x, float a, unsigned long size, cl_device_type type);
        void fill_double(void * x, double a, unsigned long size, cl_device_type type);

        template <typename DT_>
        DT_ norm_l2_false(void * x, unsigned long size, cl_device_type type, std::string function);

        template <typename DT_>
        void scale(void * x, DT_ a, unsigned long size, cl_device_type type, std::string function);

        template <typename DT_>
        void scaled_sum(void * r, void * x, void * y, DT_ b, unsigned long size, cl_device_type type, std::string function);
        void scaled_sum(void * r, void * x, void * y, unsigned long size, cl_device_type type, std::string function);

        void sum(void * r, void * x, void * y, unsigned long size, cl_device_type type, std::string function);

        template <typename DT_>
        void product_smell_dv(void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type, std::string function);

        template <typename DT_>
        void defect_smell_dv(void * rhs, void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
                unsigned long stride, unsigned long threads, cl_device_type type, std::string function);
    }
}
#endif
