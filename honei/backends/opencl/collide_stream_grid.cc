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
        void collide_stream_grid(unsigned long start, unsigned long end,
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
                DT_ tau, unsigned long size,
                cl_device_type type, std::string /*function*/)
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
            filename += "collide_stream_grid.cl";
            std::string opname("collide_stream_grid_0_");
            opname += typeid(DT_).name();
            kernel = OpenCLBackend::instance()->create_kernel(filename, opname, context, device);
            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &f_eq_0);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_0);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_temp_0);
            clSetKernelArg(kernel, 5, sizeof(DT_), (void *)&tau);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            opname = "collide_stream_grid_n_";
            opname += typeid(DT_).name();
            kernel = OpenCLBackend::instance()->create_kernel(filename, opname, context, device);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_1);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_1);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_1);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_1);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_2);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_2);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_2);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_2);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_3);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_3);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_3);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_3);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_4);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_4);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_4);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_4);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_5);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_5);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_5);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_5);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_6);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_6);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_6);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_6);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_7);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_7);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_7);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_7);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);

            clSetKernelArg(kernel, 0, sizeof(unsigned long), (void *)&start);
            clSetKernelArg(kernel, 1, sizeof(unsigned long), (void *)&end);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &dir_8);
            clSetKernelArg(kernel, 3, sizeof(cl_mem), &f_eq_8);
            clSetKernelArg(kernel, 4, sizeof(cl_mem), &f_8);
            clSetKernelArg(kernel, 5, sizeof(cl_mem), &f_temp_8);
            clSetKernelArg(kernel, 6, sizeof(DT_), (void *)&tau);
            clSetKernelArg(kernel, 7, sizeof(unsigned long), (void *)&size);
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);
        }

        template void collide_stream_grid<float> (unsigned long, unsigned long,
                void *, void *, void *, void *,
                void *, void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                float, unsigned long,
                cl_device_type, std::string);
        template void collide_stream_grid<double> (unsigned long, unsigned long,
                void *, void *, void *, void *,
                void *, void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                void *, void *, void *,
                double, unsigned long,
                cl_device_type, std::string);
    }
}
