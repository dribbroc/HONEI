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
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/backends/cuda/transfer.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>
#include <honei/util/profiler.hh>

#include <cuda_runtime.h>

using namespace honei;

namespace
{
    class cudaDefectBMQ1DVfloat
    {
        private:
            DenseVector<float> & result;
            const DenseVectorContinuousBase<float> & rhs;
            const BandedMatrixQx<Q1Type, float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaDefectBMQ1DVfloat(DenseVector<float> & result, const DenseVectorContinuousBase<float> & rhs, const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
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
                        b_gpu, result_gpu, a.size(), blocksize, a.root());

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
            }
    };

    class cudaDefectBMQ1DVdouble
    {
        private:
            DenseVector<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const BandedMatrixQx<Q1Type, double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaDefectBMQ1DVdouble(DenseVector<double> & result, const DenseVectorContinuousBase<double> & rhs, const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
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

                cuda_defect_q1_double(rhs_gpu, ll_gpu, ld_gpu, lu_gpu,
                        dl_gpu, dd_gpu, du_gpu,
                        ul_gpu, ud_gpu, uu_gpu,
                        b_gpu, result_gpu, a.size(), blocksize, a.root());

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
            }
    };

    class cudaDefectSMELLDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const DenseVectorContinuousBase<float> & rhs;
            const SparseMatrixELL<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaDefectSMELLDVfloat(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs, const SparseMatrixELL<float> & a, const DenseVectorContinuousBase<float> & b, unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_smell_dv_float(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu, b_gpu,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class cudaDefectSMELLDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const SparseMatrixELL<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaDefectSMELLDVdouble(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs, const SparseMatrixELL<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_smell_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu, b_gpu,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class MCcudaDefectSMELLDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const DenseVectorContinuousBase<float> & rhs;
            const SparseMatrixELL<float> & a;
            void * b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            MCcudaDefectSMELLDVfloat(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs, const SparseMatrixELL<float> & a, void * b, unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_smell_dv_float(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu, b,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class MCcudaDefectSMELLDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const SparseMatrixELL<double> & a;
            void * b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            MCcudaDefectSMELLDVdouble(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs, const SparseMatrixELL<double> & a, void * b, unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_smell_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu, b,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class cudaDefectSMCSRDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const DenseVectorContinuousBase<float> & rhs;
            const SparseMatrixCSR<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long atomicsize;
            unsigned long blocksize;
        public:
            cudaDefectSMCSRDVfloat(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs, const SparseMatrixCSR<float> & a, const DenseVectorContinuousBase<float> & b, unsigned long atomicsize, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                atomicsize(atomicsize),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ar_gpu(a.Ar().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_csr_dv_float(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Ar_gpu, b_gpu,
                        a.rows(), atomicsize, blocksize);

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };

    class cudaDefectSMCSRDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const SparseMatrixCSR<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long atomicsize;
            unsigned long blocksize;
        public:
            cudaDefectSMCSRDVdouble(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs, const SparseMatrixCSR<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long atomicsize, unsigned long blocksize) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                atomicsize(atomicsize),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ar_gpu(a.Ar().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_csr_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Ar_gpu, b_gpu,
                        a.rows(), atomicsize, blocksize);

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };
}

DenseVector<float> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & rhs,
        const BandedMatrixQx<Q1Type, float> & a,
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

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectBMQ1DVfloat task(result, rhs, a, b, blocksize);
        task();
    }
    else
    {
        cudaDefectBMQ1DVfloat task(result, rhs, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & rhs,
        const BandedMatrixQx<Q1Type, double> & a,
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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectBMQ1DVdouble task(result, rhs, a, b, blocksize);
        task();
    }
    else
    {
        cudaDefectBMQ1DVdouble task(result, rhs, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

DenseVector<float> & Defect<tags::GPU::CUDA>::value(DenseVector<float> & result, const DenseVectorContinuousBase<float> & rhs,
        const BandedMatrixQx<Q1Type, float> & a,
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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectBMQ1DVfloat task(result, rhs, a, b, blocksize);
        task();
    }
    else
    {
        cudaDefectBMQ1DVfloat task(result, rhs, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Defect<tags::GPU::CUDA>::value(DenseVector<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const BandedMatrixQx<Q1Type, double> & a,
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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectBMQ1DVdouble task(result, rhs, a, b, blocksize);
        task();
    }
    else
    {
        cudaDefectBMQ1DVdouble task(result, rhs, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMELLDVfloat task(result, rhs, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMELLDVfloat task(result, rhs, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> Defect<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (CUDA):");
    PROFILER_START("Defect SMELL double tags::GPU::CUDA");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMELLDVdouble task(result, rhs, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMELLDVdouble task(result, rhs, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    PROFILER_STOP("Defect SMELL double tags::GPU::CUDA");
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


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMELLDVfloat task(result, rhs, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMELLDVfloat task(result, rhs, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Defect<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixELL<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating Defect<double> (CUDA):");
    PROFILER_START("Defect SMELL double tags::GPU::CUDA");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (rhs.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(rhs.size(), a.columns());
    }


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMELLDVdouble task(result, rhs, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMELLDVdouble task(result, rhs, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    PROFILER_STOP("Defect SMELL double tags::GPU::CUDA");
    return result;
}
#endif

DenseVectorContinuousBase<float> & Defect<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs,
        const SparseMatrixCSR<float> & a,
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


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMCSRDVfloat task(result, rhs, a, b, a.blocksize(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMCSRDVfloat task(result, rhs, a, b, a.blocksize(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Defect<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
        const SparseMatrixCSR<double> & a,
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


    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDefectSMCSRDVdouble task(result, rhs, a, b, a.blocksize(), blocksize);
        task();
    }
    else
    {
        cudaDefectSMCSRDVdouble task(result, rhs, a, b, a.blocksize(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

DenseVectorContinuousBase<float> & Defect<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & rhs,
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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
        return result;
    }
    else
    {
        if (cuda::GPUPool::instance()->get_num_gpus() > 1)
        {
            //erst liegen in jeder gpu die hälfte von b vor.
            //es muss zunächst das ganze b in jede gpu gepackt werden
            //und am ende wieder das b halbiert werden | ich fasse die b haelften nur lesend an, ist also nicht nötig
            DenseVectorRange<float> b1(b.range(b.size()/2, 0));
            DenseVectorRange<float> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));

            cuda_set_device(0);
            void * whole_b_1 = cuda_malloc(b.size() * sizeof(float));
            void * whole_b_1_2 = (float*)whole_b_1 + b1.size();
            cudaMemcpy(whole_b_1, b1.lock(lm_read_only, tags::GPU::CUDA::memory_value), b1.size() * sizeof(float), cudaMemcpyDeviceToDevice);
            b1.unlock(lm_read_only);

            cuda_set_device(1);
            void * whole_b_2 = cuda_malloc(b.size() * sizeof(float));
            void * whole_b_2_2 = (float *)whole_b_2 + (b1.size());

            cudaMemcpy(whole_b_2_2 , b2.lock(lm_read_only, tags::GPU::CUDA::memory_value), b2.size() * sizeof(float), cudaMemcpyDeviceToDevice);
            b2.unlock(lm_read_only);

            cuda_set_device(0);
            cudaMemcpyPeer(whole_b_2, 1, whole_b_1, 0, b1.size() * sizeof(float));
            cudaMemcpyPeer(whole_b_1_2, 0, whole_b_2_2, 1, b2.size() * sizeof(float));

            DenseVectorRange<float> result1(result.range(result.size()/2, 0));
            DenseVectorRange<float> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
            DenseVectorRange<float> rhs1(rhs.range(rhs.size()/2, 0));
            DenseVectorRange<float> rhs2(rhs.range(rhs.size()/2 + rhs.size()%2, rhs.size()/2));
            MCcudaDefectSMELLDVfloat task1(result1, rhs1, a, whole_b_1, 0, result1.size(), blocksize);
            MCcudaDefectSMELLDVfloat task2(result2, rhs2, a, whole_b_2, result1.size(), a.rows(), blocksize);
            TicketVector tickets;
            tickets.push_back(cuda::GPUPool::instance()->enqueue(task1, 0));
            tickets.push_back(cuda::GPUPool::instance()->enqueue(task2, 1));
            tickets.wait();

            cuda_free(whole_b_1);
            cuda_set_device(1);
            cuda_free(whole_b_2);
            return result;
        }
        else
        {
            //todo DIRK remove lock hack
            b.lock(lm_read_and_write);
            b.unlock(lm_read_and_write);
            DenseVectorRange<float> result1(result.range(result.size()/2, 0));
            DenseVectorRange<float> rhs1(rhs.range(rhs.size()/2, 0));
            cudaDefectSMELLDVfloat task1(result1, rhs1, a, b, 0, result1.size(), blocksize);
            DenseVectorRange<float> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
            DenseVectorRange<float> rhs2(rhs.range(rhs.size()/2 + rhs.size()%2, rhs.size()/2));
            cudaDefectSMELLDVfloat task2(result2, rhs2, a, b, result1.size(), a.rows(), blocksize);
            cuda::GPUPool::instance()->enqueue(task1, 0).wait();
            cuda::GPUPool::instance()->enqueue(task2, 1).wait();
            b.lock(lm_read_and_write);
            b.unlock(lm_read_and_write);
            return result;
        }
    }
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Defect<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs,
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

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
        return result;
    }
    else
    {
        if (cuda::GPUPool::instance()->get_num_gpus() > 1)
        {
            //erst liegen in jeder gpu die hälfte von b vor.
            //es muss zunächst das ganze b in jede gpu gepackt werden
            //und am ende wieder das b halbiert werden | ich fasse die b haelften nur lesend an, ist also nicht nötig
            DenseVectorRange<double> b1(b.range(b.size()/2, 0));
            DenseVectorRange<double> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));

            cuda_set_device(0);
            void * whole_b_1 = cuda_malloc(b.size() * sizeof(double));
            void * whole_b_1_2 = (double*)whole_b_1 + b1.size();
            cudaMemcpy(whole_b_1, b1.lock(lm_read_only, tags::GPU::CUDA::memory_value), b1.size() * sizeof(double), cudaMemcpyDeviceToDevice);
            b1.unlock(lm_read_only);

            cuda_set_device(1);
            void * whole_b_2 = cuda_malloc(b.size() * sizeof(double));
            void * whole_b_2_2 = (double *)whole_b_2 + (b1.size());

            cudaMemcpy(whole_b_2_2 , b2.lock(lm_read_only, tags::GPU::CUDA::memory_value), b2.size() * sizeof(double), cudaMemcpyDeviceToDevice);
            b2.unlock(lm_read_only);

            cuda_set_device(0);
            cudaMemcpyPeer(whole_b_2, 1, whole_b_1, 0, b1.size() * sizeof(double));
            cudaMemcpyPeer(whole_b_1_2, 0, whole_b_2_2, 1, b2.size() * sizeof(double));

            DenseVectorRange<double> result1(result.range(result.size()/2, 0));
            DenseVectorRange<double> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
            DenseVectorRange<double> rhs1(rhs.range(rhs.size()/2, 0));
            DenseVectorRange<double> rhs2(rhs.range(rhs.size()/2 + rhs.size()%2, rhs.size()/2));
            MCcudaDefectSMELLDVdouble task1(result1, rhs1, a, whole_b_1, 0, result1.size(), blocksize);
            MCcudaDefectSMELLDVdouble task2(result2, rhs2, a, whole_b_2, result1.size(), a.rows(), blocksize);
            TicketVector tickets;
            tickets.push_back(cuda::GPUPool::instance()->enqueue(task1, 0));
            tickets.push_back(cuda::GPUPool::instance()->enqueue(task2, 1));
            tickets.wait();

            cuda_free(whole_b_1);
            cuda_set_device(1);
            cuda_free(whole_b_2);
            return result;
        }
        else
        {
            //todo DIRK remove lock hack
            b.lock(lm_read_and_write);
            b.unlock(lm_read_and_write);
            DenseVectorRange<double> result1(result.range(result.size()/2, 0));
            DenseVectorRange<double> rhs1(rhs.range(rhs.size()/2, 0));
            cudaDefectSMELLDVdouble task1(result1, rhs1, a, b, 0, result1.size(), blocksize);
            DenseVectorRange<double> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
            DenseVectorRange<double> rhs2(rhs.range(rhs.size()/2 + rhs.size()%2, rhs.size()/2));
            cudaDefectSMELLDVdouble task2(result2, rhs2, a, b, result1.size(), a.rows(), blocksize);
            cuda::GPUPool::instance()->enqueue(task1, 0).wait();
            cuda::GPUPool::instance()->enqueue(task2, 1).wait();
            b.lock(lm_read_and_write);
            b.unlock(lm_read_and_write);
            return result;
        }
    }
}
#endif
