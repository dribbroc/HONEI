/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
        DT_ dot_product(void * x, void * y, unsigned long size, cl_device_type type, std::string function)
        {
            cl_command_queue command_queue;
            cl_kernel kernel;
            cl_context context;
            cl_device_id device;

            DCQ dcq = OpenCLBackend::instance()->prepare_device(type);
            device = dcq.device;
            context = dcq.context;
            command_queue = dcq.command_queue;


            //print_device_info(device);
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/backends/opencl/";
            filename += "dot_product.cl";
            kernel = OpenCLBackend::instance()->create_kernel(filename, function, context, device);
            size_t tmp_work_group_size;
            clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void*)&tmp_work_group_size, NULL);
            size_t blocksize = tmp_work_group_size;
            size_t threads = blocksize * 32;
            cl_mem tmp_device(0);
            DT_ * tmp_cpu = new DT_[threads];
            tmp_device = OpenCLBackend::instance()->create_empty_buffer(threads*sizeof(DT_), dcq.context);
            clSetKernelArg(kernel, 0, sizeof(cl_mem), &x);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), &y);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), &tmp_device);
            clSetKernelArg(kernel, 3, sizeof(cl_uint), (void *)&size);

            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, &blocksize, 0, NULL, NULL);
            clEnqueueReadBuffer(dcq.command_queue, tmp_device, CL_TRUE, 0, threads*sizeof(DT_), (void*)tmp_cpu, 0, NULL, NULL);
            clReleaseMemObject(tmp_device);
            clFinish(command_queue);
            DT_ result(0);
            for (unsigned long i(0) ; i < threads ; ++i)
            {
                result += tmp_cpu[i];
            }
            delete[] tmp_cpu;

            return result;
        }
        template float dot_product<float>(void*, void*, unsigned long, cl_device_type, std::string);
        template double dot_product<double>(void*, void*, unsigned long, cl_device_type, std::string);
    }
}
