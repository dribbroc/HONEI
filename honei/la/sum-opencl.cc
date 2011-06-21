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

#include <honei/la/sum.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>
#include <typeinfo>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_sum(DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
        {
            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            void * x_cl(x.lock(lm_read_and_write, Tag_::memory_value));
            void * y_cl(y.lock(lm_read_only, Tag_::memory_value));
            std::string opname("sum_three_");
            opname += typeid(DT_).name();
            sum(x_cl, x_cl, y_cl, x.size(), tag_to_device<Tag_>(), opname);
            x.unlock(lm_read_and_write);
            y.unlock(lm_read_only);
        }
    }
}

using namespace honei;

template <typename DT_>
DenseVectorContinuousBase<DT_> & Sum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating the sum of two DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("Sum tags::OpenCL::CPU");
    opencl::common_sum<tags::OpenCL::CPU>(x, y);
    PROFILER_STOP("Sum tags::OpenCL::CPU");
    return x;
}
template DenseVectorContinuousBase<float> & Sum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template DenseVectorContinuousBase<double> & Sum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);

template <typename DT_>
DenseVectorContinuousBase<DT_> & Sum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating the sum of two DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("Sum tags::OpenCL::GPU");
    opencl::common_sum<tags::OpenCL::GPU>(x, y);
    PROFILER_STOP("Sum tags::OpenCL::GPU");
    return x;
}
template DenseVectorContinuousBase<float> & Sum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template DenseVectorContinuousBase<double> & Sum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);
