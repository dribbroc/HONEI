/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2013 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

namespace honei
{
    namespace opencl
    {
        template <typename DT_>
        void extraction_grid_wet(unsigned long start, unsigned long end,
                void * f_0, void * f_1, void * f_2,
                void * f_3, void * f_4, void * f_5,
                void * f_6, void * f_7, void * f_8,
                void * h, void * u, void * v,
                void * distribution_x, void * distribution_y, DT_ epsilon,
                cl_device_type type, std::string function)
        {
            cl_command_queue command_queue;
            cl_kernel kernel;
            cl_context context;
            cl_device_id device;
            size_t threads = end - start;

            DCQ dcq = OpenCLBackend::instance()->prepare_device(type);
            device = dcq.device;
            context = dcq.context;
            command_queue = dcq.command_queue;

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/backends/opencl/";
            filename += "extraction_grid.cl";
            kernel = OpenCLBackend::instance()->create_kernel(filename, function, context, device);
            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &f_0);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_1);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_2);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_3);
            clSetKernelArg(kernel, 6, sizeof(cl_mem), &f_4);
            clSetKernelArg(kernel, 7, sizeof(cl_mem), &f_5);
            clSetKernelArg(kernel, 8, sizeof(cl_mem), &f_6);
            clSetKernelArg(kernel, 9, sizeof(cl_mem), &f_7);
            clSetKernelArg(kernel, 10, sizeof(cl_mem), &f_8);
            clSetKernelArg(kernel, 11, sizeof(cl_mem), &h);
            clSetKernelArg(kernel, 12, sizeof(cl_mem), &u);
            clSetKernelArg(kernel, 13, sizeof(cl_mem), &v);
            clSetKernelArg(kernel, 14, sizeof(cl_mem), &distribution_x);
            clSetKernelArg(kernel, 15, sizeof(cl_mem), &distribution_y);
            clSetKernelArg(kernel, 16, sizeof(DT_), (void *)&epsilon);
            clSetKernelArg(kernel, 17, 18*sizeof(DT_), NULL);

            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);
        }

        template void extraction_grid_wet<float> (unsigned long, unsigned long,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, float,
                cl_device_type, std::string);
        template void extraction_grid_wet<double> (unsigned long, unsigned long,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, double,
                cl_device_type, std::string);
    }
}
