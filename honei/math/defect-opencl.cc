/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008, 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/defect.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

using namespace honei;

DenseVector<float> Defect<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (OpenCL CPU):");
    PROFILER_START("Defect tags::OpenCL::CPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::defect_smell_dv_float(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::CPU");
    return result;
}

DenseVector<double> Defect<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (OpenCL CPU):");
    PROFILER_START("Defect tags::OpenCL::CPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::defect_smell_dv_double(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::CPU");
    return result;
}

DenseVector<float> Defect<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (OpenCL GPU):");
    PROFILER_START("Defect tags::OpenCL::GPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::defect_smell_dv_float(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::GPU");
    return result;
}

DenseVector<double> Defect<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (OpenCL GPU):");
    PROFILER_START("Defect tags::OpenCL::GPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::defect_smell_dv_double(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::GPU");
    return result;
}

DenseVectorContinuousBase<float> & Defect<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (OpenCL CPU):");
    PROFILER_START("Defect tags::OpenCL::CPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::defect_smell_dv_float(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::CPU");
    return result;
}

DenseVectorContinuousBase<double> & Defect<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (OpenCL CPU):");
    PROFILER_START("Defect tags::OpenCL::CPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::defect_smell_dv_double(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::CPU");
    return result;
}

DenseVectorContinuousBase<float> & Defect<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (OpenCL GPU):");
    PROFILER_START("Defect tags::OpenCL::GPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::defect_smell_dv_float(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::GPU");
    return result;
}

DenseVectorContinuousBase<double> & Defect<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (OpenCL GPU):");
    PROFILER_START("Defect tags::OpenCL::GPU");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    void * rhs_cl(rhs.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::defect_smell_dv_double(rhs_cl, b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    PROFILER_STOP("Defect tags::OpenCL::GPU");
    return result;
}
