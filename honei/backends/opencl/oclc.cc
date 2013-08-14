/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of HONEI. HONEI is free software;
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


#include <iostream>
#include <honei/util/file_to_string.hh>
#include <honei/util/stringify.hh>
#include <honei/backends/opencl/opencl_backend.hh>

using namespace honei;
using namespace honei::opencl;

void compile(cl_device_type type, std::string filename)
{
    DCQ dcq = OpenCLBackend::instance()->prepare_device(type);
    cl_device_id device = dcq.device;
    cl_context context = dcq.context;

    std::string file(HONEI_SOURCEDIR);
    file+= "/honei/backends/opencl/";
    file+= filename;
    //kernel = OpenCLBackend::instance()->create_kernel(filename, function, context, device);

    std::string  sourceStr = file_to_string(file);
    const char * source = sourceStr.c_str();
    size_t source_size[] = { strlen(source) };
    cl_program program(NULL);
    cl_int status(CL_SUCCESS);

    program = clCreateProgramWithSource(context, 1, &source, source_size, &status);

    if (status != CL_SUCCESS)
    {
        throw InternalError("OpenCL: Error " + stringify(status) + " in clCreateProgramWithSource!");
    }
    std::string cl_flags("-Werror -cl-strict-aliasing -cl-mad-enable -cl-fast-relaxed-math");

    if (type == CL_DEVICE_TYPE_CPU)
        cl_flags += " -DHONEI_OPENCL_CPU";
    else if (type == CL_DEVICE_TYPE_GPU)
        cl_flags += " -DHONEI_CUDA_DOUBLE -DHONEI_OPENCL_GPU";
    else if (type == CL_DEVICE_TYPE_ACCELERATOR)
        cl_flags += " -DHONEI_OPENCL_ACCELERATOR";

    status = clBuildProgram(program, 1, &device, cl_flags.c_str(), NULL, NULL);
    if (status != CL_SUCCESS)
    {
        OpenCLBackend::instance()->print_device_info(device);
        OpenCLBackend::instance()->print_program_info(program, device);
        throw InternalError("OpenCL: Error " + stringify(status) + " in clBuildProgram!");
    }
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cout<<"Usage 'ellinfo ocl file'"<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string filename(argv[1]);
    compile(CL_DEVICE_TYPE_CPU , filename);
    compile(CL_DEVICE_TYPE_GPU , filename);

    std::cout<<"OCLC compiled "<<argv[1]<<" successfully"<<std::endl;
    return EXIT_SUCCESS;
}
