/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#include <honei/la/scaled_sum.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_scaled_sum(DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y, DT_ b)
        {
            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            void * x_cl(x.lock(lm_read_and_write, Tag_::memory_value));
            void * y_cl(y.lock(lm_read_only, Tag_::memory_value));
            std::string opname("scaled_sum_three_s_");
            opname += typeid(DT_).name();
            scaled_sum(x_cl, x_cl, y_cl, b, x.size(), tag_to_device<Tag_>(), opname);
            x.unlock(lm_read_and_write);
            y.unlock(lm_read_only);
        }

        template <typename Tag_, typename DT_>
        void common_scaled_sum(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y, DT_ b)
        {
            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());
            if (r.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            void * r_cl(r.lock(lm_write_only, Tag_::memory_value));
            void * x_cl(x.lock(lm_read_only, Tag_::memory_value));
            void * y_cl(y.lock(lm_read_only, Tag_::memory_value));
            std::string opname("scaled_sum_three_s_");
            opname += typeid(DT_).name();
            scaled_sum(r_cl, x_cl, y_cl, b, x.size(), tag_to_device<Tag_>(), opname);
            r.unlock(lm_write_only);
            x.unlock(lm_read_only);
            y.unlock(lm_read_only);
        }

        template <typename Tag_, typename DT_>
        void common_scaled_sum(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
        {
            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());
            if (r.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            void * r_cl(r.lock(lm_read_and_write, Tag_::memory_value));
            void * x_cl(x.lock(lm_read_only, Tag_::memory_value));
            void * y_cl(y.lock(lm_read_only, Tag_::memory_value));
            std::string opname("scaled_sum_three_");
            opname += typeid(DT_).name();
            opencl::scaled_sum(r_cl, x_cl, y_cl, x.size(), tag_to_device<Tag_>(), opname);
            r.unlock(lm_read_and_write);
            x.unlock(lm_read_only);
            y.unlock(lm_read_only);
        }
    }
}

using namespace honei;

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");
    opencl::common_scaled_sum<tags::OpenCL::CPU>(x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return x;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");
    opencl::common_scaled_sum<tags::OpenCL::GPU>(x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return x;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL Accelerator):");
    PROFILER_START("ScaledSum tags::OpenCL::Accelerator");
    opencl::common_scaled_sum<tags::OpenCL::Accelerator>(x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::Accelerator");
    return x;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");
    opencl::common_scaled_sum<tags::OpenCL::CPU>(r, x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");
    opencl::common_scaled_sum<tags::OpenCL::GPU>(r, x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y, DT_ b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL Accelerator):");
    PROFILER_START("ScaledSum tags::OpenCL::Accelerator");
    opencl::common_scaled_sum<tags::OpenCL::Accelerator>(r, x, y, b);
    PROFILER_STOP("ScaledSum tags::OpenCL::Accelerator");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");
    opencl::common_scaled_sum<tags::OpenCL::CPU>(r, x, y);
    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");
    opencl::common_scaled_sum<tags::OpenCL::GPU>(r, x, y);
    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);

template <typename DT_>
DenseVectorContinuousBase<DT_> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x,
        const DenseVectorContinuousBase<DT_> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<DT_> (OpenCL Accelerator):");
    PROFILER_START("ScaledSum tags::OpenCL::Accelerator");
    opencl::common_scaled_sum<tags::OpenCL::Accelerator>(r, x, y);
    PROFILER_STOP("ScaledSum tags::OpenCL::Accelerator");
    return r;
}
template DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &);
template DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::Accelerator>::value(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &);
