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

#include <honei/la/product.hh>
#include <honei/backends/opencl/operations.hh>

using namespace honei;

DenseVector<float> Product<tags::OpenCL::CPU>::value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<float> with DenseVectorContinuousBase<float> (OpenCL CPU):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::product_smell_dv_float(b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    return result;
}

DenseVector<double> Product<tags::OpenCL::CPU>::value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<double> with DenseVectorContinuousBase<double> (OpenCL CPU):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    void * b_cl(b.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::product_smell_dv_double(b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), CL_DEVICE_TYPE_CPU);
    result.unlock(lm_write_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    return result;
}

DenseVector<float> Product<tags::OpenCL::GPU>::value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<float> with DenseVectorContinuousBase<float> (OpenCL GPU):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::product_smell_dv_float(b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    return result;
}

DenseVector<double> Product<tags::OpenCL::GPU>::value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<double> with DenseVectorContinuousBase<double> (OpenCL GPU):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    void * b_cl(b.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * result_cl(result.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * Aj_cl(a.Aj().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Ax_cl(a.Ax().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * Arl_cl(a.Arl().lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::product_smell_dv_double(b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), CL_DEVICE_TYPE_GPU);
    result.unlock(lm_write_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);
    a.Arl().unlock(lm_read_only);

    return result;
}
