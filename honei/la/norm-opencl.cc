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

#include <honei/la/norm.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
            DT_ common_norm_l2_false(const DenseVectorContinuousBase<DT_> & x)
            {
                DT_ result(0);
                if (x.size() < 32 * 128)
                {
                    x.lock(lm_read_only);
                    for (unsigned long i(0) ; i < x.size() ; ++i)
                    {
                        result += x[i] * x[i];
                    }
                    x.unlock(lm_read_only);
                }
                else
                {
                    void * x_cl(x.lock(lm_read_only, Tag_::memory_value));
                    std::string opname("norm_l2_false_");
                    opname += typeid(DT_).name();
                    result = norm_l2_false<DT_>(x_cl, x.size(), tag_to_device<Tag_>(), opname);
                    x.unlock(lm_read_only);
                }
                return result;
            }
    }
}

using namespace honei;

template <typename DT_>
DT_ Norm<vnt_l_two, false, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::CPU");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::CPU>(x);
    PROFILER_STOP("Norm l2 false tags::OpenCL::CPU");
    return result;
}
template float Norm<vnt_l_two, false, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, false, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ Norm<vnt_l_two, true, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::CPU");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::CPU>(x);
    PROFILER_STOP("Norm l2 true tags::OpenCL::CPU");
    return sqrt(result);
}
template float Norm<vnt_l_two, true, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, true, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ Norm<vnt_l_two, false, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::GPU");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::GPU>(x);
    PROFILER_STOP("Norm l2 false tags::OpenCL::GPU");
    return result;
}
template float Norm<vnt_l_two, false, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, false, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ Norm<vnt_l_two, true, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::GPU");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::GPU>(x);
    PROFILER_STOP("Norm l2 true tags::OpenCL::GPU");
    return sqrt(result);
}
template float Norm<vnt_l_two, true, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, true, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ Norm<vnt_l_two, false, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL Accelerator):");
    PROFILER_START("Norm l2 false tags::OpenCL::Accelerator");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::Accelerator>(x);
    PROFILER_STOP("Norm l2 false tags::OpenCL::Accelerator");
    return result;
}
template float Norm<vnt_l_two, false, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, false, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ Norm<vnt_l_two, true, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<DT_> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<DT_> (OpenCL Accelerator):");
    PROFILER_START("Norm l2 true tags::OpenCL::Accelerator");
    DT_ result = opencl::common_norm_l2_false<tags::OpenCL::Accelerator>(x);
    PROFILER_STOP("Norm l2 true tags::OpenCL::Accelerator");
    return sqrt(result);
}
template float Norm<vnt_l_two, true, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<float> &);
template double Norm<vnt_l_two, true, tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<double> &);
