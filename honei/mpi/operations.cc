/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/tags.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>
#include <honei/mpi/sparse_matrix_csr_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/norm.hh>
#include <honei/la/product.hh>
#include <honei/math/defect.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>
#include <honei/la/element_product.hh>

#ifdef HONEI_CUDA
#include <honei/backends/cuda/transfer.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <honei/util/time_stamp.hh>

namespace
{
#ifdef HONEI_CUDA
    template <typename DT_>
    class AllocTask
    {
        private:
            void ** address;
            unsigned long size;
        public:
            AllocTask(void** a, unsigned long s) :
                address(a),
                size(s)
        {
        }

            void operator() ()
            {
                *address = cuda_malloc_host(size);
            }
    };

    template <typename DT_>
    class DownloadTask
    {
        private:
            const DenseVector<DT_> & device;
            void * address;
            cudaStream_t stream;
        public:
            DownloadTask(const DenseVector<DT_> & d, void * a, cudaStream_t s = 0) :
                device(d),
                address(a),
                stream(s)
        {
        }

            void operator() ()
            {
                void * d_gpu(device.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_download_async(d_gpu, address, device.size() * sizeof(DT_), stream);
                /// \todo this is a realy bad idea!
                cudaStreamSynchronize(stream);
                device.unlock(lm_read_only);
            }
    };

    template <typename DT_>
    class UploadTask
    {
        private:
            const DenseVector<DT_> & device;
            void * address;
            cudaStream_t stream;
        public:
            UploadTask(void * a, const DenseVector<DT_> & d, cudaStream_t s = 0) :
                device(d),
                address(a),
                stream(s)
        {
        }

            void operator() ()
            {
                void * d_gpu(device.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                cuda_upload_async(address, d_gpu, device.size() * sizeof(DT_), stream);
                device.unlock(lm_write_only);
            }
    };

    template <typename MT_>
    class ProductTask
    {
    };

    template <>
    class ProductTask<SparseMatrixELLMPI<double> >
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const SparseMatrixELL<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            ProductTask(DenseVectorContinuousBase<double> & result, const SparseMatrixELL<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize, cudaStream_t s) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize),
                stream(s)
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
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads(), stream);

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    template <>
    class ProductTask<SparseMatrixCSRMPI<double> >
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const SparseMatrixCSR<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long atomicsize;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            ProductTask(DenseVectorContinuousBase<double> & result, const SparseMatrixCSR<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize, cudaStream_t s) :
                result(result),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize),
                stream(s)
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
                        row_start, row_end, a.blocksize(), blocksize, stream);

                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };

    template <typename MT_>
    class DefectTask
    {
    };

    template <>
    class DefectTask<SparseMatrixELLMPI<double> >
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const SparseMatrixELL<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            DefectTask(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs, const SparseMatrixELL<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize, cudaStream_t s) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize),
                stream(s)
            {
            }

            void operator() ()
            {
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * rhs_gpu(rhs.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * Aj_gpu(a.Aj().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Ax_gpu(a.Ax().lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * Arl_gpu(a.Arl().lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_defect_smell_dv_double(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, Arl_gpu, b_gpu,
                        row_start, row_end, a.num_cols_per_row(), a.stride(), blocksize, a.threads(), stream);

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Arl().unlock(lm_read_only);
            }
    };

    template <>
    class DefectTask<SparseMatrixCSRMPI<double> >
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & rhs;
            const SparseMatrixCSR<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long row_start;
            unsigned long row_end;
            unsigned long atomicsize;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            DefectTask(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & rhs, const SparseMatrixCSR<double> & a, const DenseVectorContinuousBase<double> & b,
                    unsigned long row_start, unsigned long row_end, unsigned long blocksize, cudaStream_t s) :
                result(result),
                rhs(rhs),
                a(a),
                b(b),
                row_start(row_start),
                row_end(row_end),
                blocksize(blocksize),
                stream(s)
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
                        a.rows(), a.blocksize(), blocksize, stream);

                result.unlock(lm_write_only);
                rhs.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.Aj().unlock(lm_read_only);
                a.Ax().unlock(lm_read_only);
                a.Ar().unlock(lm_read_only);
            }
    };

    class SumTask
    {
        private:
            DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            SumTask(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize, cudaStream_t s) :
                a(a),
                b(b),
                blocksize(blocksize),
                stream(s)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_sum_two_double(a_gpu, b_gpu, a.size(), blocksize, stream);
                cudaDeviceSynchronize();

                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class DifferenceTask
    {
        private:
            DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
            cudaStream_t stream;
        public:
            DifferenceTask(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize, cudaStream_t s) :
                a(a),
                b(b),
                blocksize(blocksize),
                stream(s)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_two_double(a_gpu, b_gpu, a.size(), blocksize, stream);
                cudaDeviceSynchronize();

                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };
#endif
}

using namespace honei;

static void * temp_data = 0;
static unsigned long temp_data_size = 0;

static void * temp_send_data = 0;
static unsigned long temp_send_data_size = 0;

static void * temp_missing_data = 0;
static unsigned long temp_missing_data_size = 0;

#ifdef HONEI_CUDA
static cudaStream_t * streams = 0;
#endif

template <typename Tag_>
    template <typename DT_>
void MPIOps<Tag_>::difference(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Difference<Tag_>::value(r.vector(), x.vector(), y.vector());
}

#ifdef HONEI_CUDA
    template <typename DT_>
void MPIOps<tags::GPU::CUDA>::difference(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Difference<tags::GPU::CUDA>::value(r.vector(), x.vector(), y.vector());
}
#endif

template <typename Tag_>
    template <typename DT_>
DT_ MPIOps<Tag_>::dot_product(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    DT_ local_result(DotProduct<Tag_>::value(x.vector(), y.vector()));
    DT_ result(DT_(0));

    MPI_Allreduce(&local_result, &result, 1, mpi::MPIType<DT_>::value(), MPI_SUM, MPI_COMM_WORLD);

    return result;
}

#ifdef HONEI_CUDA
    template <typename DT_>
DT_ MPIOps<tags::GPU::CUDA>::dot_product(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    DT_ local_result(DotProduct<tags::GPU::CUDA>::value(x.vector(), y.vector()));
    DT_ result(DT_(0));

    MPI_Allreduce(&local_result, &result, 1, mpi::MPIType<DT_>::value(), MPI_SUM, MPI_COMM_WORLD);

    return result;
}
#endif

template <typename Tag_>
    template <typename DT_>
void MPIOps<Tag_>::element_product(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    ElementProduct<Tag_>::value(r.vector(), x.vector(), y.vector());
}

#ifdef HONEI_CUDA
    template <typename DT_>
void MPIOps<tags::GPU::CUDA>::element_product(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    ElementProduct<tags::GPU::CUDA>::value(r.vector(), x.vector(), y.vector());
}
#endif

template <typename Tag_>
    template <typename DT_>
DT_ MPIOps<Tag_>::norm_l2_false(const DenseVectorMPI<DT_> & x)
{
    DT_ local_result(Norm<vnt_l_two, false, Tag_>::value(x.vector()));
    DT_ result(DT_(0));

    MPI_Allreduce(&local_result, &result, 1, mpi::MPIType<DT_>::value(), MPI_SUM, MPI_COMM_WORLD);

    return result;
}

#ifdef HONEI_CUDA
    template <typename DT_>
DT_ MPIOps<tags::GPU::CUDA>::norm_l2_false(const DenseVectorMPI<DT_> & x)
{
    DT_ local_result(Norm<vnt_l_two, false, tags::GPU::CUDA>::value(x.vector()));
    DT_ result(DT_(0));

    MPI_Allreduce(&local_result, &result, 1, mpi::MPIType<DT_>::value(), MPI_SUM, MPI_COMM_WORLD);

    return result;
}
#endif

template <typename Tag_>
template <typename MT_, typename DT_>
void MPIOps<Tag_>::product(DenseVectorMPI<DT_> & r, const MT_ & a, const DenseVectorMPI<DT_> & b)
{

    // berechne innere anteile
    TicketVector inner_ticket;
    if (a.active())
    {
        OperationWrapper<honei::Product<Tag_>, DenseVector<DT_>,
            DenseVector<DT_>, typename MT_::LocalType_, DenseVector<DT_>, unsigned long, unsigned long > wrapper(r.vector());
        inner_ticket.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, r.vector(), a.inner_matrix(), b.vector(), 0, a.inner_matrix().rows())));
    }

    const DT_ * bp = b.elements();

    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    // empfange alle fehlenden werte
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values.elements() + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    if (temp_send_data_size < a.send_size() * sizeof(DT_))
    {
        ::free(temp_send_data);
        temp_send_data = ::malloc(a.send_size() * sizeof(DT_));
        temp_send_data_size = a.send_size() * sizeof(DT_);
    }

    DT_ * send_data((DT_*)temp_send_data);
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = bp[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    // berechne innere anteile
    //if (a.active()) Product<Tag_>::value(r.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    inner_ticket.wait();
    if (a.active()) Product<Tag_>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Sum<Tag_>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
}

#ifdef HONEI_CUDA
template <typename MT_, typename DT_>
void MPIOps<tags::GPU::CUDA>::product(DenseVectorMPI<DT_> & r, const MT_ & a, const DenseVectorMPI<DT_> & b)
{
    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    unsigned long blocksize_prod(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));
    unsigned long blocksize_sum(Configuration::instance()->get_value("cuda::sum_two_double", 128ul));
    TicketVector tickets_0;
    TicketVector tickets_1;

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    if (streams == NULL)
    {
        streams = new cudaStream_t[2];
        //inner product stream
        cudaStreamCreate(&streams[0]);
        // outer/synchronisation stream
        cudaStreamCreate(&streams[1]);
    }

    // berechne innere anteile
    if (a.active())
    {
        ProductTask<MT_> product_inner(r.vector(), a.inner_matrix(), b.vector(), 0, a.inner_matrix().rows(), blocksize_prod, streams[0]);
        tickets_0.push_back(cuda::GPUPool::instance()->enqueue(product_inner, 0));
    }

    // empfange alle fehlenden werte
    if (temp_missing_data_size < a.outer_matrix().columns() * sizeof(DT_))
    {
        cuda_free_host(temp_missing_data);
        AllocTask<DT_> at(&temp_missing_data, a.outer_matrix().columns() * sizeof(DT_));
        cuda::GPUPool::instance()->enqueue(at, 0).wait();
        temp_missing_data_size = a.outer_matrix().columns() * sizeof(DT_);
    }
    DT_ * missing_values_array((DT_*)temp_missing_data);
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values_array + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    if (temp_send_data_size < a.send_size() * sizeof(DT_))
    {
        ::free(temp_send_data);
        temp_send_data = ::malloc(a.send_size() * sizeof(DT_));
        temp_send_data_size = a.send_size() * sizeof(DT_);
    }
    DT_ * send_data((DT_*)temp_send_data);

    if (temp_data_size < b.local_size() * sizeof(DT_))
    {
        cuda_free_host(temp_data);
        AllocTask<DT_> at(&temp_data, b.local_size() * sizeof(DT_));
        cuda::GPUPool::instance()->enqueue(at, 0).wait();
        temp_data_size = b.local_size() * sizeof(DT_);
    }

    DownloadTask<DT_> dt(b.vector(), temp_data, streams[1]);
    cuda::GPUPool::instance()->enqueue(dt, 0).wait();
    DT_ * b_cpu((DT_*) temp_data);


    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = b_cpu[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    //fehlende aeussere werte auf gpu schaufeln
    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    UploadTask<DT_> ut(temp_missing_data, missing_values, streams[1]);
    tickets_1.push_back(cuda::GPUPool::instance()->enqueue(ut, 0));

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active())
    {
        ProductTask<MT_> product_outer(r_outer, a.outer_matrix(), missing_values, 0, a.outer_matrix().rows(), blocksize_prod, streams[1]);
        tickets_1.push_back(cuda::GPUPool::instance()->enqueue(product_outer, 0));
        SumTask sum_outer(r.vector(), r_outer, blocksize_sum, streams[1]);
        tickets_1.push_back(cuda::GPUPool::instance()->enqueue(sum_outer, 0));
    }

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    tickets_0.wait();
    tickets_1.wait();
    //cudaStreamSynchronize(streams[0]);
    //cudaStreamSynchronize(streams[1]);
    //cudaDeviceSynchronize();
}
#endif

template <typename Tag_>
template <typename MT_, typename DT_>
void MPIOps<Tag_>::defect(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & rhs, const MT_ & a, const DenseVectorMPI<DT_> & b)
{
    // berechne innere anteile
    TicketVector inner_ticket;
    if (a.active())
    {
        OperationWrapper<honei::Defect<Tag_>, DenseVector<DT_>,
            DenseVector<DT_>, DenseVector<DT_>, typename MT_::LocalType_, DenseVector<DT_>, unsigned long, unsigned long > wrapper(r.vector());
        inner_ticket.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, r.vector(), rhs.vector(), a.inner_matrix(), b.vector(), 0, a.inner_matrix().rows())));
    }

    const DT_ * bp = b.elements();

    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    // empfange alle fehlenden werte
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values.elements() + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    if (temp_send_data_size < a.send_size() * sizeof(DT_))
    {
        ::free(temp_send_data);
        temp_send_data = ::malloc(a.send_size() * sizeof(DT_));
        temp_send_data_size = a.send_size() * sizeof(DT_);
    }
    DT_ * send_data((DT_*)temp_send_data);
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = bp[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    // berechne innere anteile
    //if (a.active()) Defect<Tag_>::value(r.vector(), rhs.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    inner_ticket.wait();
    if (a.active()) Product<Tag_>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Difference<Tag_>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
}

#ifdef HONEI_CUDA
template <typename MT_, typename DT_>
void MPIOps<tags::GPU::CUDA>::defect(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & rhs, const MT_ & a, const DenseVectorMPI<DT_> & b)
{
    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    unsigned long blocksize_prod(Configuration::instance()->get_value("cuda::product_smell_dv_double", 128ul));
    unsigned long blocksize_sum(Configuration::instance()->get_value("cuda::sum_two_double", 128ul));
    TicketVector tickets_0;
    TicketVector tickets_1;

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    if (streams == NULL)
    {
        streams = new cudaStream_t[2];
        //inner product stream
        cudaStreamCreate(&streams[0]);
        // outer/synchronisation stream
        cudaStreamCreate(&streams[1]);
    }

    // berechne innere anteile
    if (a.active())
    {
        DefectTask<MT_> defect_inner(r.vector(), rhs.vector(), a.inner_matrix(), b.vector(), 0, a.inner_matrix().rows(), blocksize_prod, streams[0]);
        tickets_0.push_back(cuda::GPUPool::instance()->enqueue(defect_inner, 0));
    }

    // empfange alle fehlenden werte
    if (temp_missing_data_size < a.outer_matrix().columns() * sizeof(DT_))
    {
        cuda_free_host(temp_missing_data);
        AllocTask<DT_> at(&temp_missing_data, a.outer_matrix().columns() * sizeof(DT_));
        cuda::GPUPool::instance()->enqueue(at, 0).wait();
        temp_missing_data_size = a.outer_matrix().columns() * sizeof(DT_);
    }
    DT_ * missing_values_array((DT_*)temp_missing_data);
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values_array + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    if (temp_send_data_size < a.send_size() * sizeof(DT_))
    {
        ::free(temp_send_data);
        temp_send_data = ::malloc(a.send_size() * sizeof(DT_));
        temp_send_data_size = a.send_size() * sizeof(DT_);
    }
    DT_ * send_data((DT_*)temp_send_data);

    if (temp_data_size < b.local_size() * sizeof(DT_))
    {
        cuda_free_host(temp_data);
        AllocTask<DT_> at(&temp_data, b.local_size() * sizeof(DT_));
        cuda::GPUPool::instance()->enqueue(at, 0).wait();
        temp_data_size = b.local_size() * sizeof(DT_);
    }

    DownloadTask<DT_> dt(b.vector(), temp_data, streams[1]);
    cuda::GPUPool::instance()->enqueue(dt, 0).wait();
    DT_ * b_cpu((DT_*) temp_data);


    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = b_cpu[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    //fehlende aeussere werte auf gpu schaufeln
    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    UploadTask<DT_> ut(temp_missing_data, missing_values, streams[1]);
    tickets_1.push_back(cuda::GPUPool::instance()->enqueue(ut, 0));

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active())
    {
        ProductTask<MT_> product_outer(r_outer, a.outer_matrix(), missing_values, 0, a.outer_matrix().rows(), blocksize_prod, streams[1]);
        tickets_1.push_back(cuda::GPUPool::instance()->enqueue(product_outer, 0));
        DifferenceTask difference_outer(r.vector(), r_outer, blocksize_sum, streams[1]);
        tickets_1.push_back(cuda::GPUPool::instance()->enqueue(difference_outer, 0));
    }

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    tickets_0.wait();
    tickets_1.wait();
    //cudaStreamSynchronize(streams[0]);
    //cudaStreamSynchronize(streams[1]);
    //cudaDeviceSynchronize();
}
#endif

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::scale(DenseVectorMPI<DT_> & x, DT_ a)
{
    Scale<Tag_>::value(x.vector(), a);
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::scale(DenseVectorMPI<DT_> & x, DT_ a)
{
    Scale<tags::GPU::CUDA>::value(x.vector(), a);
}
#endif

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::scaled_sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a)
{
    ScaledSum<Tag_>::value(x.vector(), y.vector(), a);
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a)
{
    ScaledSum<tags::GPU::CUDA>::value(x.vector(), y.vector(), a);
}
#endif

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::scaled_sum(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a)
{
    ScaledSum<Tag_>::value(r.vector(), x.vector(), y.vector(), a);
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a)
{
    ScaledSum<tags::GPU::CUDA>::value(r.vector(), x.vector(), y.vector(), a);
}
#endif

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Sum<Tag_>::value(x.vector(), y.vector());
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Sum<tags::GPU::CUDA>::value(x.vector(), y.vector());
}
#endif

template <typename Tag_>
template <typename MT_>
MT_ MPIOps<Tag_>::transposition(const MT_ & src)
{
    typedef typename MT_::DataType DT_;

    src.inner_matrix().lock(lm_read_only);
    src.inner_matrix().unlock(lm_read_only);
    src.outer_matrix().lock(lm_read_only);
    src.outer_matrix().unlock(lm_read_only);

    // reverse create src local matrix
    unsigned long new_rows(DenseVectorMPI<double>::calc_size(src.columns()));
    SparseMatrix<DT_> local(src.local_rows(), src.columns());
    SparseMatrix<DT_> old_inner(src.inner_matrix());
    for(typename SparseMatrix<DT_>::NonZeroConstElementIterator i(old_inner.begin_non_zero_elements()) ; i != old_inner.end_non_zero_elements() ; ++i)
    {
        local(i.row(), i.column() + src.x_offset(), *i);
    }
    SparseMatrix<DT_> old_outer(src.outer_matrix());
    unsigned long cix(0);
    for (std::set<unsigned long>::iterator ci(src.missing_indices().begin()) ; ci != src.missing_indices().end() ; ++ci, ++cix)
    {
        for (unsigned long i(0) ; i < old_outer.column(cix).used_elements() ; ++i)
        {
            local(old_outer.column(cix).indices()[i], *ci, old_outer.column(cix).elements()[i]);
        }
    }

    // create new transposed local matrix
    SparseMatrix<DT_> new_local(new_rows, src.rows());
    unsigned long start_row(DenseVectorMPI<double>::calc_offset(src.columns()));
    unsigned long end_row(start_row + new_rows);

    unsigned long _rank(mpi::mpi_comm_rank());
    unsigned long _com_size(mpi::mpi_comm_size());
    // jeder schickt nach einander alle elemente einer spalte->zeile an den aktuellen empfänger
    for (unsigned long rank(0) ; rank < _com_size ; ++rank)
    {
        if (rank == _rank)
        {
            mpi::mpi_bcast(&start_row, 1, rank);
            mpi::mpi_bcast(&end_row, 1, rank);

            for (unsigned long row(start_row) ; row < end_row ; ++row)
            {
                for (unsigned long other(0) ; other < _com_size ; ++other)
                {
                    if (other == _rank)
                    {
                        for (unsigned long i(0) ; i < local.column(row).used_elements() ; ++i)
                        {
                            new_local(row - start_row, local.column(row).indices()[i] + src.offset(), local.column(row).elements()[i]);
                        }
                    }
                    else
                    {
                        unsigned long count(0);
                        mpi::mpi_recv(&count, 1, other, other);

                        if (count > 0)
                        {
                            unsigned long alien_offset;
                            mpi::mpi_recv(&alien_offset, 1, other, other);

                            unsigned long col_indices[count];
                            DT_ values[count];

                            mpi::mpi_recv(col_indices, count, other, other);
                            mpi::mpi_recv(values, count, other, other);

                            for (unsigned long i(0) ; i < count ; ++i)
                            {
                                new_local(row - start_row, col_indices[i] + alien_offset, values[i]);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            unsigned long alien_start_row;
            unsigned long alien_end_row;
            mpi::mpi_bcast(&alien_start_row, 1, rank);
            mpi::mpi_bcast(&alien_end_row, 1, rank);

            for (unsigned long col(alien_start_row) ; col < alien_end_row ; ++col)
            {
                unsigned long count(local.column(col).used_elements());
                mpi::mpi_send(&count, 1, rank, _rank);

                if (count > 0)
                {
                    unsigned long own_offset(src.offset());
                    mpi::mpi_send(&own_offset, 1, rank, _rank);
                    mpi::mpi_send(local.column(col).indices(), count, rank, _rank);
                    mpi::mpi_send(local.column(col).elements(), count, rank, _rank);
                }
            }
        }
    }

    MT_ result(new_local, src.columns());
    return result;
}


template struct MPIOps<tags::CPU>;
template void MPIOps<tags::CPU>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template SparseMatrixELLMPI<double> MPIOps<tags::CPU>::transposition(const SparseMatrixELLMPI<double> & src);
template SparseMatrixCSRMPI<double> MPIOps<tags::CPU>::transposition(const SparseMatrixCSRMPI<double> & src);

template struct MPIOps<tags::CPU::Generic>;
template void MPIOps<tags::CPU::Generic>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::Generic>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::Generic>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::Generic>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::Generic>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::Generic>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::Generic>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::Generic>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::Generic>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::Generic>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::Generic>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::Generic>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

template struct MPIOps<tags::CPU::MultiCore::Generic>;
template void MPIOps<tags::CPU::MultiCore::Generic>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::Generic>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::Generic>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::MultiCore::Generic>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::Generic>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::MultiCore::Generic>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::MultiCore::Generic>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::Generic>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::Generic>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::MultiCore::Generic>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::Generic>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::Generic>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

#ifdef HONEI_SSE
template struct MPIOps<tags::CPU::SSE>;
template void MPIOps<tags::CPU::SSE>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::SSE>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::SSE>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::SSE>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::SSE>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::SSE>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::SSE>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::SSE>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

template struct MPIOps<tags::CPU::MultiCore::SSE>;
template void MPIOps<tags::CPU::MultiCore::SSE>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::MultiCore::SSE>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::SSE>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::MultiCore::SSE>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::MultiCore::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::SSE>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
#endif

#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
template struct MPIOps<tags::GPU::CUDA>;
template void MPIOps<tags::GPU::CUDA>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::GPU::CUDA>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::GPU::CUDA>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::GPU::CUDA>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::GPU::CUDA>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::GPU::CUDA>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::GPU::CUDA>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::GPU::CUDA>::product(DenseVectorMPI<double> & r, const SparseMatrixCSRMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::GPU::CUDA>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::GPU::CUDA>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
#endif
#endif
