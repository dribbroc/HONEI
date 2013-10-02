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

#include <honei/la/dot_product.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
            DT_ common_dot_product(const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                DT_ result(0);
                if (x.size() < 32 * 128)
                {
                    x.lock(lm_read_only);
                    y.lock(lm_read_only);
                    for (unsigned long i(0) ; i < x.size() ; ++i)
                    {
                        result += x[i] * y[i];
                    }
                    x.unlock(lm_read_only);
                    y.unlock(lm_read_only);
                }
                else
                {
                    void * x_cl(x.lock(lm_read_only, Tag_::memory_value));
                    void * y_cl(y.lock(lm_read_only, Tag_::memory_value));
                    std::string opname("dot_product_");
                    opname += typeid(DT_).name();
                    result = dot_product<DT_>(x_cl, y_cl, x.size(), tag_to_device<Tag_>(), opname);
                    x.unlock(lm_read_only);
                    y.unlock(lm_read_only);
                }
                return result;
            }
    }
}

using namespace honei;

template <typename DT_>
DT_ DotProduct<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<DT_>(OpenCL CPU):");
    PROFILER_START("DotProduct tags::OpenCL::CPU");
    DT_ result = opencl::common_dot_product<tags::OpenCL::CPU>(x, y);
    PROFILER_STOP("DotProduct tags::OpenCL::CPU");
    return result;
}
template float DotProduct<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template double DotProduct<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ DotProduct<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<DT_>(OpenCL GPU):");
    PROFILER_START("DotProduct tags::OpenCL::GPU");
    DT_ result = opencl::common_dot_product<tags::OpenCL::GPU>(x, y);
    PROFILER_STOP("DotProduct tags::OpenCL::GPU");
    return result;
}
template float DotProduct<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template double DotProduct<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);

template <typename DT_>
DT_ DotProduct<tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<DT_>(OpenCL Accelerator):");
    PROFILER_START("DotProduct tags::OpenCL::Accelerator");
    DT_ result = opencl::common_dot_product<tags::OpenCL::Accelerator>(x, y);
    PROFILER_STOP("DotProduct tags::OpenCL::Accelerator");
    return result;
}
template float DotProduct<tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template double DotProduct<tags::OpenCL::Accelerator>::value(const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);
