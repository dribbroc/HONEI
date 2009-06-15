/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>

using namespace honei;

DenseVector<float> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & rhs,
        const BandedMatrixQ1<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_float", 128ul));
    unsigned long m(a.root());

    void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * ll_gpu(a.band(LL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ld_gpu(a.band(LD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * lu_gpu(a.band(LU).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * dl_gpu(a.band(DL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * dd_gpu(a.band(DD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * du_gpu(a.band(DU).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ul_gpu(a.band(UL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ud_gpu(a.band(UD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * uu_gpu(a.band(UU).lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_defect_q1_float(rhs_gpu, ll_gpu, ld_gpu, lu_gpu,
            dl_gpu, dd_gpu, du_gpu,
            ul_gpu, ud_gpu, uu_gpu,
            b_gpu, result_gpu, a.size(), blocksize, m);

    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.band(LL).unlock(lm_read_only);
    a.band(LD).unlock(lm_read_only);
    a.band(LU).unlock(lm_read_only);
    a.band(DL).unlock(lm_read_only);
    a.band(DD).unlock(lm_read_only);
    a.band(DU).unlock(lm_read_only);
    a.band(UL).unlock(lm_read_only);
    a.band(UD).unlock(lm_read_only);
    a.band(UU).unlock(lm_read_only);

    return result;
}

DenseVector<float> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 256ul));

    void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_defect_smell_dv_float(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), blocksize);

    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 256ul));

    void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_defect_smell_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), blocksize);

    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);

    return result;
}
#endif

DenseVectorContinuousBase<float> & Defect<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixELL<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating Defect<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 256ul));

    void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_defect_smell_dv_float(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), blocksize);

    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Defect<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 256ul));

    void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_defect_smell_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), blocksize);

    result.unlock(lm_write_only);
    rhs.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.Aj().unlock(lm_read_only);
    a.Ax().unlock(lm_read_only);

    return result;
}
#endif
