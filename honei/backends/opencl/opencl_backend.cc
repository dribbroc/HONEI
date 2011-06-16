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

#include <honei/backends/opencl/opencl_backend.hh>
#include <honei/backends/opencl/opencl_backend-impl.hh>
#include <honei/util/instantiation_policy-impl.hh>

namespace honei
{
    template class InstantiationPolicy<OpenCLBackend, Singleton>;
    namespace opencl
    {
        class OpenCLBackend;

        template <typename Tag_>
        cl_device_type tag_to_device()
        {
            throw InternalError("Device not supported!");
            return 0;
        }

        template <>
        cl_device_type tag_to_device<tags::OpenCL::CPU>()
        {
            return CL_DEVICE_TYPE_CPU;
        }

        template <>
        cl_device_type tag_to_device<tags::OpenCL::GPU>()
        {
            return CL_DEVICE_TYPE_GPU;
        }

        template <>
        cl_device_type tag_to_device<tags::OpenCL::Accelerator>()
        {
            return CL_DEVICE_TYPE_ACCELERATOR;
        }
    }
}
