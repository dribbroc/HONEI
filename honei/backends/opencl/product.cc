/* vim: set sw=4 sts=4 et nofoldenable : */

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

#include <honei/backends/opencl/opencl_backend.hh>
#include <honei/util/attributes.hh>
#include <iostream>

namespace honei
{
    namespace opencl
    {
        template <typename DT_>
        void product_smell_dv(void * x, void * y, void * Aj, void * Ax, void * Arl,
                unsigned long num_rows, HONEI_UNUSED unsigned long num_cols, HONEI_UNUSED unsigned long num_cols_per_row,
                unsigned long stride, unsigned long ell_threads, cl_device_type type, std::string function)
        {
            cl_command_queue command_queue;
            cl_kernel kernel;
            cl_context context;
            cl_device_id device;
            size_t threads;
            if (type == CL_DEVICE_TYPE_CPU)
                threads = num_rows;
            else
                threads = num_rows * ell_threads;

            DCQ dcq = OpenCLBackend::instance()->prepare_device(type);
            device = dcq.device;
            context = dcq.context;
            command_queue = dcq.command_queue;

            //print_device_info(device);
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/backends/opencl/";
            filename += "product.cl";
            kernel = OpenCLBackend::instance()->create_kernel(filename, function, context, device);
            clSetKernelArg(kernel, 0, sizeof(cl_mem), &x);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), &y);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &Aj);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &Ax);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &Arl);
            clSetKernelArg(kernel, 5, sizeof(unsigned long), (void *)&num_rows);
            clSetKernelArg(kernel, 6, sizeof(unsigned long), (void *)&stride);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&ell_threads);
            size_t tmp_work_group_size;
            clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void*)&tmp_work_group_size, NULL);
            clSetKernelArg(kernel, 8, tmp_work_group_size*sizeof(DT_), NULL);

            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);
        }
        template void product_smell_dv<float>(void*, void*, void*, void*, void*,
                unsigned long, unsigned long, unsigned long, unsigned long, unsigned long, cl_device_type, std::string);
        template void product_smell_dv<double>(void*, void*, void*, void*, void*,
                unsigned long, unsigned long, unsigned long, unsigned long, unsigned long, cl_device_type, std::string);
    }
}
