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

#include <honei/util/unittest.hh>
#include <honei/backends/opencl/opencl_backend.hh>
#include <iostream>
#include <honei/backends/opencl/operations.hh>

using namespace honei;
using namespace honei::opencl;
using namespace tests;

template <typename Tag_>
class OpenclPlatformQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        OpenclPlatformQuickTest() :
            QuickTaggedTest<Tag_>("opencl_platform_test")
        {
        }

        virtual void run() const
        {
            OpenCLBackend::instance()->print_platform_info();
        }
};
OpenclPlatformQuickTest<tags::OpenCL> opencl_platform_quick_test;

template <typename Tag_>
class OpenclBackendQuickTest :
    public QuickTaggedTest<Tag_>
{
    private:
        cl_device_type _type;
    public:
        OpenclBackendQuickTest(cl_device_type type) :
            QuickTaggedTest<Tag_>("opencl_backend_test")
        {
            _type = type;
        }

        virtual void run() const
        {
            cl_command_queue command_queue(NULL);
            cl_context context(NULL);
            cl_device_id device(NULL);
            cl_kernel kernel(NULL);
            cl_mem x_buffer;
            cl_mem y_buffer;

            DCQ dcq = OpenCLBackend::instance()->prepare_device(_type);
            device = dcq.device;
            context = dcq.context;
            command_queue = dcq.command_queue;
            OpenCLBackend::instance()->print_device_info(device);
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/backends/opencl/";
            filename += "scaled_sum.cl";
            kernel = OpenCLBackend::instance()->create_kernel(filename, "scaled_sum_three_s_f", context, device);

            float x[10];
            float y[10];
            unsigned long size = 10;
            for (unsigned long i(0) ; i < size ; ++i)
            {
                x[i] = i;
                y[i] = 2;
            }
            for (unsigned long i(0) ; i < size ; ++i)
                std::cout<<x[i]<<" ";
            std::cout<<std::endl;
            //x_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(cl_float) * size, x, NULL);
            //y_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cl_float) * size, y, NULL);
            x_buffer = OpenCLBackend::instance()->create_buffer(sizeof(cl_float) * size, context, x);
            y_buffer = OpenCLBackend::instance()->create_buffer(sizeof(cl_float) * size, context, y);

            float b = 2;
            clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&x_buffer);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&x_buffer);
            clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&y_buffer);
            clSetKernelArg(kernel, 3, sizeof(cl_float), (void *)&b);
            clSetKernelArg(kernel, 4, sizeof(unsigned long), (void *)&size);

            size_t threads = size;
            clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &threads, NULL, 0, NULL, NULL);
            clFinish(command_queue);

            clEnqueueReadBuffer(command_queue, x_buffer, CL_TRUE, 0, size * sizeof(cl_float), x, 0, NULL, NULL);
            clFinish(command_queue);
            for (unsigned long i(0) ; i < size ; ++i)
            {
                std::cout<<x[i]<<" ";
                TEST_CHECK_EQUAL(x[i], float(i + 4));
            }
            std::cout<<std::endl;
        }
};
OpenclBackendQuickTest<tags::OpenCL::CPU> cpu_opencl_backend_quick_test(CL_DEVICE_TYPE_CPU);
#ifdef HONEI_OPENCL_GPU
OpenclBackendQuickTest<tags::OpenCL::GPU> gpu_opencl_backend_quick_test(CL_DEVICE_TYPE_GPU);
#endif
#ifdef HONEI_OPENCL_ACC
OpenclBackendQuickTest<tags::OpenCL::Accelerator> accelerator_opencl_backend_quick_test(CL_DEVICE_TYPE_ACCELERATOR);
#endif
