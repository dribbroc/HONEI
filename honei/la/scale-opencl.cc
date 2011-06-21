/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/scale.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_scale(DenseVectorContinuousBase<DT_> & x, const DT_ a)
        {
            void * x_cl(x.lock(lm_read_and_write, Tag_::memory_value));
            std::string opname("scale_");
            opname += typeid(DT_).name();
            scale(x_cl, a, x.size(), tag_to_device<Tag_>(), opname);
            x.unlock(lm_read_and_write);
        }
    }
}

using namespace honei;

template <typename DT_>
DenseVectorContinuousBase<DT_> & Scale<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<DT_> & x, const DT_ a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<DT_> by scalar (OpenCL CPU):");
    PROFILER_START("Scale tags::OpenCL::CPU");
    opencl::common_scale<tags::OpenCL::CPU>(x, a);
    PROFILER_STOP("Scale tags::OpenCL::CPU");
    return x;
}
template DenseVectorContinuousBase<float> & Scale<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> &, const float);
template DenseVectorContinuousBase<double> & Scale<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> &, const double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & Scale<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<DT_> & x, const DT_ a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<DT_> by scalar (OpenCL GPU):");
    PROFILER_START("Scale tags::OpenCL::GPU");
    opencl::common_scale<tags::OpenCL::GPU>(x, a);
    PROFILER_STOP("Scale tags::OpenCL::GPU");
    return x;
}
template DenseVectorContinuousBase<float> & Scale<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> &, const float);
template DenseVectorContinuousBase<double> & Scale<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> &, const double);
