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

#include <honei/la/product.hh>
#include <honei/la/algorithm.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>
#include <honei/util/profiler.hh>

using namespace honei;

namespace
{
    class cudaProductBMDVfloat
    {
        private:
            void * result_gpu;
            void * band_gpu;
            void * temp_b_gpu;
            unsigned long size;
            unsigned long blocksize;
        public:
            cudaProductBMDVfloat(void * result_gpu, void * band_gpu, void * temp_b_gpu, unsigned long size, unsigned long blocksize) :
                result_gpu(result_gpu),
                band_gpu(band_gpu),
                temp_b_gpu(temp_b_gpu),
                size(size),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                cuda_scaled_sum_three_float(result_gpu, band_gpu, temp_b_gpu, size, blocksize);
            }
    };

    class cudaProductBMQ1DVfloat
    {
        private:
            DenseVector<float> & result;;
            const BandedMatrixQx<Q1Type, float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaProductBMQ1DVfloat(DenseVector<float> & result, const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
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

                cuda_product_bmdv_q1_float(ll_gpu, ld_gpu, lu_gpu,
                        dl_gpu, dd_gpu, du_gpu,
                        ul_gpu, ud_gpu, uu_gpu,
                        b_gpu, result_gpu, a.size(), blocksize, a.root());

                result.unlock(lm_write_only);
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

    class cudaProductBMQ1DVdouble
    {
        private:
            DenseVector<double> & result;;
            const BandedMatrixQx<Q1Type, double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaProductBMQ1DVdouble(DenseVector<double> & result, const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
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

                cuda_product_bmdv_q1_double(ll_gpu, ld_gpu, lu_gpu,
                        dl_gpu, dd_gpu, du_gpu,
                        ul_gpu, ud_gpu, uu_gpu,
                        b_gpu, result_gpu, a.size(), blocksize, a.root());

                result.unlock(lm_write_only);
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

    class cudaProductSMELLDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const SparseMatrixELL<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaProductSMELLDVfloat(DenseVectorContinuousBase<float> & result, const SparseMatrixELL<float> & a, const DenseVectorContinuousBase<float> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_product_smell_dv_float(b_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class cudaProductSMELLDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const SparseMatrixELL<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaProductSMELLDVdouble(DenseVectorContinuousBase<double> & result, const SparseMatrixELL<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_product_smell_dv_double(b_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads());

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    class cudaProductSMCSRDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const SparseMatrixCSR<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaProductSMCSRDVfloat(DenseVectorContinuousBase<float> & result, const SparseMatrixCSR<float> & a, const DenseVectorContinuousBase<float> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ar_gpu(a.Ar().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_product_csr_dv_float(b_gpu, result_gpu, Aj_gpu, Ax_gpu, Ar_gpu,
                        row_start, row_end, blocksize);

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };

    class cudaProductSMCSRDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const SparseMatrixCSR<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
        public:
            cudaProductSMCSRDVdouble(DenseVectorContinuousBase<double> & result, const SparseMatrixCSR<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ar_gpu(a.Ar().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_product_csr_dv_double(b_gpu, result_gpu, Aj_gpu, Ax_gpu, Ar_gpu,
                        row_start, row_end, blocksize);

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };
}

DenseVector<float> Product<tags::GPU::CUDA>::value(const BandedMatrix<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows(), float(0));

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_float", 128ul));

    unsigned long middle_index(a.rows() - 1);
    unsigned long op_offset;

    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

    for (BandedMatrix<float>::ConstBandIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            void * band_gpu(band->lock(lm_read_only, tags::GPU::CUDA::memory_value));
            op_offset = band.index() - middle_index;
            void * temp_b_gpu = (void *)((float *)b_gpu + op_offset);
            if (! cuda::GPUPool::instance()->idle())
            {
                cudaProductBMDVfloat task(result_gpu, band_gpu, temp_b_gpu, a.size() - op_offset, blocksize);
                task();
            }
            else
            {
                cudaProductBMDVfloat task(result_gpu, band_gpu, temp_b_gpu, a.size() - op_offset, blocksize);
                cuda::GPUPool::instance()->enqueue(task, 0).wait();
            }
            band->unlock(lm_read_only);
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            void * band_gpu(band->lock(lm_read_only, tags::GPU::CUDA::memory_value));
            band_gpu = (void *)((float *)band_gpu + op_offset);
            void * temp_result_gpu = (void *)((float *)result_gpu + op_offset);
            if (! cuda::GPUPool::instance()->idle())
            {
                cudaProductBMDVfloat task(temp_result_gpu, band_gpu, b_gpu, a.size() - op_offset, blocksize);
                task();
            }
            else
            {
                cudaProductBMDVfloat task(temp_result_gpu, band_gpu, b_gpu, a.size() - op_offset, blocksize);
                cuda::GPUPool::instance()->enqueue(task, 0).wait();
            }
            band->unlock(lm_read_only);
        }
    }

    result.unlock(lm_read_and_write);
    b.unlock(lm_read_only);

    return result;
}


DenseVector<float> Product<tags::GPU::CUDA>::value(const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrixQx<Q1Type, float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductBMQ1DVfloat task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaProductBMQ1DVfloat task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> Product<tags::GPU::CUDA>::value(const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying BandedMatrixQx<Q1Type, double> with DenseVectorContinuousBase<double> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductBMQ1DVdouble task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaProductBMQ1DVdouble task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

DenseVector<float> & Product<tags::GPU::CUDA>::value(DenseVector<float> & result, const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrixQx<Q1Type, float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (result.size() != a.rows())
    {
        throw VectorSizeDoesNotMatch(result.size(), a.rows());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductBMQ1DVfloat task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaProductBMQ1DVfloat task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Product<tags::GPU::CUDA>::value(DenseVector<double> & result, const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying BandedMatrixQx<Q1Type, double> with DenseVectorContinuousBase<double> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (result.size() != a.rows())
    {
        throw VectorSizeDoesNotMatch(result.size(), a.rows());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductBMQ1DVdouble task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaProductBMQ1DVdouble task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

DenseVector<float> & Product<tags::GPU::CUDA>::value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductSMELLDVfloat task(result, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaProductSMELLDVfloat task(result, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Product<tags::GPU::CUDA>::value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<double> with DenseVectorContinuousBase<double> (CUDA):");
    PROFILER_START("Product SMELL double tags::GPU::CUDA");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductSMELLDVdouble task(result, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaProductSMELLDVdouble task(result, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    PROFILER_STOP("Product SMELL double tags::GPU::CUDA");
    return result;
}
#endif

DenseVector<float> & Product<tags::GPU::CUDA>::value(DenseVector<float> & result, const SparseMatrixCSR<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying SparseMatrixCSR<float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductSMCSRDVfloat task(result, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaProductSMCSRDVfloat task(result, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Product<tags::GPU::CUDA>::value(DenseVector<double> & result, const SparseMatrixCSR<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying SparseMatrixCSR<double> with DenseVectorContinuousBase<double> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProductSMCSRDVdouble task(result, a, b, 0, a.rows(), blocksize);
        task();
    }
    else
    {
        cudaProductSMCSRDVdouble task(result, a, b, 0, a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return result;
}
#endif

DenseVector<float> & Product<tags::GPU::MultiCore::CUDA>::value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<float> with DenseVectorContinuousBase<float> (MC CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        //todo DIRK remove lock hack
        b.lock(lm_read_and_write);
        b.unlock(lm_read_and_write);
        DenseVectorRange<float> result1(result.range(result.size()/2, 0));
        cudaProductSMELLDVfloat task1(result1, a, b, 0, result1.size(), blocksize);
        DenseVectorRange<float> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        cudaProductSMELLDVfloat task2(result2, a, b, result1.size(), a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0).wait();
        cuda::GPUPool::instance()->enqueue(task2, 1).wait();
        b.lock(lm_read_and_write);
        b.unlock(lm_read_and_write);
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Product<tags::GPU::MultiCore::CUDA>::value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<double> with DenseVectorContinuousBase<double> (MC CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }
    if (a.rows() != result.size())
    {
        throw VectorSizeDoesNotMatch(a.rows(), result.size());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        //todo DIRK remove lock hack
        b.lock(lm_read_and_write);
        b.unlock(lm_read_and_write);
        DenseVectorRange<double> result1(result.range(result.size()/2, 0));
        cudaProductSMELLDVdouble task1(result1, a, b, 0, result1.size(), blocksize);
        DenseVectorRange<double> result2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        cudaProductSMELLDVdouble task2(result2, a, b, result1.size(), a.rows(), blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0).wait();
        cuda::GPUPool::instance()->enqueue(task2, 1).wait();
        b.lock(lm_read_and_write);
        b.unlock(lm_read_and_write);
    }

    return result;
}
#endif
