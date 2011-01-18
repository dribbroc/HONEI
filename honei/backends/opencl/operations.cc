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

namespace honei
{
    namespace opencl
    {
        /*void copy_float(void * x, void * y, unsigned long size, cl_device_type type)
        {
            cl_command_queue command_queue;
            cl_kernel kernel;
            cl_context context;
            cl_device_id device;
            size_t threads =size;

            DCQ dcq = OpenCLBackend::instance()->prepare_device(type);
            device = dcq.device;
            context = dcq.context;
            command_queue = dcq.command_queue;

            //print_device_info(device);
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/backends/opencl/";
            filename += "operations.cl";
            kernel = OpenCLBackend::instance()->create_kernel(filename, "copy_float", context, device);
            clSetKernelArg(kernel, 0, sizeof(cl_mem), &x);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), &y);
            clSetKernelArg(kernel, 2, sizeof(cl_uint), (void *)&size);

            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);
            clFinish(command_queue);
        }*/
    }
}
