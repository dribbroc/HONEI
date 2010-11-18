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
#ifndef MEMORY_GUARD_OPENCL_BACKEND_IMPL_HH
#define MEMORY_GUARD_OPENCL_BACKEND_IMPL_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/backends/opencl/opencl_backend.hh>
#include <honei/util/file_to_string.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <vector>

namespace honei
{
    template <> struct Implementation<OpenCLBackend>
    {
        struct KSD
        {
            cl_kernel kernel;
            std::string kernel_name;
            cl_device_id device;
        };

        std::vector<DCQ> dcqs;
        std::vector<KSD> kernels;

        /// Our mutex.
        Mutex * const _mutex;

        /// Constructor
        Implementation():
            _mutex(new Mutex)
        {
        }

        /// Destructor
        ~Implementation()
        {
            delete _mutex;
        }

        void print_program_info(cl_program program, cl_device_id device)
        {
            Lock l(*_mutex);
            char output2[10000];
            size_t size = 0;
            clGetProgramInfo(program, CL_PROGRAM_SOURCE, 0, NULL, &size);
            char  output[size];
            clGetProgramInfo(program, CL_PROGRAM_SOURCE, size, (void *) output, NULL);
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 10000, (void *) output2, NULL);
            std::cout<<"Program Source: "<<output<<std::endl;
            std::cout<<"Build Log: "<<output2<<std::endl;
        }

        void print_device_info(cl_device_id device)
        {
            Lock l(*_mutex);
            char output[1000];
            char output2[1000];
            char output3[10000];
            clGetDeviceInfo(device, CL_DEVICE_VENDOR, 1000, (void *) output, NULL);
            std::cout<<"Device vendor: "<<output<<std::endl;
            clGetDeviceInfo(device, CL_DEVICE_NAME, 1000, (void *) output2, NULL);
            std::cout<<"Device name: "<<output2<<std::endl;
            clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 10000, (void *) output3, NULL);
            std::cout<<"Extensions: "<<output3<<std::endl;
        }

        void print_platform_info()
        {
            cl_uint numPlatforms;
            cl_platform_id platform = NULL;
            clGetPlatformIDs(0, NULL, &numPlatforms);
            if(numPlatforms > 0)
            {
                cl_platform_id* platforms = new cl_platform_id[numPlatforms];
                clGetPlatformIDs(numPlatforms, platforms, NULL);
                for(unsigned int i=0; i < numPlatforms; ++i)
                {
                    if (i != 0) std::cout<<"--------------------"<<std::endl;

                    char pbuff[1000];
                    clGetPlatformInfo(
                            platforms[i],
                            CL_PLATFORM_VENDOR,
                            sizeof(pbuff),
                            pbuff,
                            NULL);
                    platform = platforms[i];
                    std::cout<<"Platform:"<<std::endl;
                    std::cout<<pbuff<<std::endl;
                    clGetPlatformInfo(
                            platforms[i],
                            CL_PLATFORM_NAME,
                            sizeof(pbuff),
                            pbuff,
                            NULL);
                    std::cout<<pbuff<<std::endl;
                    cl_device_id devices [10];
                    cl_uint device_count(0);
                    clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 10, devices, &device_count);
                    std::cout<<"Devices:"<<std::endl;
                    for (unsigned long i(0) ; i < device_count ; ++i)
                    {
                        print_device_info(devices[i]);
                    }
                }
                delete platforms;
            }
            else throw InternalError("OpenCL: No platform found!");
        }


        cl_kernel add_kernel(std::string file, std::string kernel_name, cl_context context, cl_device_id device)
        {
            std::string  sourceStr = file_to_string(file);
            const char * source = sourceStr.c_str();
            size_t source_size[] = { strlen(source) };
            cl_program program(NULL);
            cl_int status(CL_SUCCESS);
            program = clCreateProgramWithSource(context, 1, &source, source_size, &status);
            if (status != CL_SUCCESS)
            {
                throw InternalError("OpenCL: Error " + stringify(status) + " in clCrateProgramWithSource!");
            }
            std::string cl_flags("-Werror");
            status = clBuildProgram(program, 1, &device, cl_flags.c_str(), NULL, NULL);
            if (status != CL_SUCCESS)
            {
                print_device_info(device);
                print_program_info(program, device);
                throw InternalError("OpenCL: Error " + stringify(status) + " in clBuildProgram!");
            }
            cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &status);
            if (status != CL_SUCCESS)
                throw InternalError("OpenCL: Error " + stringify(status) + " in clCreateKernel!");
            KSD temp;
            temp.kernel = kernel;
            temp.kernel_name = kernel_name;
            temp.device = device;
            kernels.push_back(temp);
            return kernel;
        }

        cl_kernel create_kernel(std::string file, std::string kernel_name, cl_context context, cl_device_id device)
        {
            Lock l(*_mutex);
            for (std::vector<KSD>::iterator i(kernels.begin()) ; i < kernels.end() ; ++i)
            {
                if (i->kernel_name == kernel_name && i->device == device)
                {
                    return i->kernel;
                }
            }
            return add_kernel(file, kernel_name, context, device);
        }

        DCQ add_device(cl_device_type type)
        {
            cl_device_id device;
            cl_context context;
            cl_command_queue command_queue;
            cl_uint numPlatforms;
            cl_platform_id platform = NULL;
            device = NULL;
            clGetPlatformIDs(0, NULL, &numPlatforms);
            if(numPlatforms > 0)
            {
                cl_platform_id* platforms = new cl_platform_id[numPlatforms];
                clGetPlatformIDs(numPlatforms, platforms, NULL);
                for(unsigned int i=0; i < numPlatforms; ++i)
                {
                    char pbuff[100];
                    clGetPlatformInfo(
                            platforms[i],
                            CL_PLATFORM_VENDOR,
                            sizeof(pbuff),
                            pbuff,
                            NULL);
                    platform = platforms[i];
                    cl_device_id devices [10];
                    cl_uint device_count(0);
                    clGetDeviceIDs(platform, type, 10, devices, &device_count);
                    if (device_count > 0)
                    {
                        device = devices[0];
                        break;
                    }
                }
                delete platforms;
            }
            else throw InternalError("OpenCL: No platform found!");
            if (device == NULL)
                throw InternalError("OpenCL: No device found!");

            cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };
            cl_context_properties* cprops = (NULL == platform) ? NULL : cps;

            context = clCreateContext(cprops, 1, &device, NULL, NULL, NULL);

            command_queue = clCreateCommandQueue( context, device, 0, NULL);
            DCQ temp;
            temp.device = device;
            temp.context = context;
            temp.command_queue = command_queue;
            temp.type = type;
            dcqs.push_back(temp);
            return temp;
        }

        DCQ prepare_device(cl_device_type type)
        {
            Lock l(*_mutex);
            for (std::vector<DCQ>::iterator i(dcqs.begin()) ; i < dcqs.end() ; ++i)
            {
                if (i->type == type)
                {
                    return *i;
                }
            }
            DCQ temp = add_device(type);
            return temp;
        }

        cl_mem create_empty_buffer(unsigned long bytes, cl_context context)
        {
            Lock l(*_mutex);
            return clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, NULL);
        }
        cl_mem create_buffer(unsigned long bytes, cl_context context, void * src)
        {
            Lock l(*_mutex);
            return clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, bytes, src, NULL);
        }
    };

    OpenCLBackend::OpenCLBackend() :
        PrivateImplementationPattern<OpenCLBackend, Shared>(new Implementation<OpenCLBackend>())
    {
    }

    OpenCLBackend::~OpenCLBackend()
    {
    }

    cl_kernel OpenCLBackend::create_kernel(std::string file, std::string kernel_name, cl_context context, cl_device_id device)
    {
        return _imp->create_kernel(file, kernel_name, context, device);
    }

    DCQ OpenCLBackend::prepare_device(cl_device_type type)
    {
        return _imp->prepare_device(type);
    }

    void OpenCLBackend::print_device_info(cl_device_id device)
    {
        _imp->print_device_info(device);
    }

    void OpenCLBackend::print_platform_info()
    {
        _imp->print_platform_info();
    }

    void OpenCLBackend::print_program_info(cl_program program, cl_device_id device)
    {
        _imp->print_program_info(program, device);
    }

    cl_mem OpenCLBackend::create_empty_buffer(unsigned long bytes, cl_context context)
    {
        return _imp->create_empty_buffer(bytes, context);
    }

    cl_mem OpenCLBackend::create_buffer(unsigned long bytes, cl_context context, void * src)
    {
        return _imp->create_buffer(bytes, context, src);
    }
}
#endif
