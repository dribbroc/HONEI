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
#endif

using namespace honei;

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
template <typename DT_>
void MPIOps<Tag_>::product(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
{
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
    DT_ * send_data = new DT_[a.send_size()];
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = bp[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    // berechne innere anteile
    if (a.active()) Product<Tag_>::value(r.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active()) Product<Tag_>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Sum<Tag_>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    delete[] send_data;
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::product(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
{
    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    DT_ * missing_values_array = (DT_*)cuda_malloc_host(a.outer_matrix().columns() * sizeof(DT_));
    // empfange alle fehlenden werte
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values_array + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    DT_ * send_data = new DT_[a.send_size()];
    DT_ * b_cpu = (DT_*)cuda_malloc_host(b.local_size() * sizeof(DT_));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_download(b_gpu, b_cpu, b.local_size() * sizeof(DT_));
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = b_cpu[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }
    b.unlock(lm_read_only);

    // berechne innere anteile
    if (a.active()) Product<tags::GPU::CUDA>::value(r.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();
    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    void * missing_values_gpu(missing_values.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    cuda_upload(missing_values_array, missing_values_gpu, missing_values.size() * sizeof(DT_));
    missing_values.unlock(lm_write_only);

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active()) Product<tags::GPU::CUDA>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Sum<tags::GPU::CUDA>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    cuda_free_host(missing_values_array);
    cuda_free_host(b_cpu);
    delete[] send_data;
}
#endif

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::defect(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & rhs, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
{
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
    DT_ * send_data = new DT_[a.send_size()];
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = bp[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }

    // berechne innere anteile
    if (a.active()) Defect<Tag_>::value(r.vector(), rhs.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active()) Product<Tag_>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Difference<Tag_>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    delete[] send_data;
}

#ifdef HONEI_CUDA
template <typename DT_>
void MPIOps<tags::GPU::CUDA>::defect(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & rhs, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
{
    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    DT_ * missing_values_array = (DT_*)cuda_malloc_host(a.outer_matrix().columns() * sizeof(DT_));
    // empfange alle fehlenden werte
    unsigned long g_size(0);
    for (unsigned long i(0) ; i < a.recv_ranks().size() ; ++i)
    {
        recv_requests.push_back(mpi::mpi_irecv(missing_values_array + g_size, a.recv_sizes().at(i), a.recv_ranks().at(i), a.recv_ranks().at(i)));
        g_size += a.recv_sizes().at(i);
    }

    // sende alle werte, die anderen fehlen
    g_size = 0;
    DT_ * send_data = new DT_[a.send_size()];
    DT_ * b_cpu = (DT_*)cuda_malloc_host(b.local_size() * sizeof(DT_));
    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_download(b_gpu, b_cpu, b.local_size() * sizeof(DT_));
    for (unsigned long i(0) ; i < a.send_ranks().size() ; ++i)
    {
        unsigned long g_end(g_size + a.send_sizes().at(i));
        for (unsigned long j(0) ; g_size < g_end ; ++g_size, ++j)
            send_data[g_size] = b_cpu[a.send_index().at(g_size)];
        send_requests.push_back(mpi::mpi_isend(&(send_data[g_size - a.send_sizes().at(i)]), a.send_sizes().at(i), a.send_ranks().at(i), myrank));
    }
    b.unlock(lm_read_only);

    // berechne innere anteile
    if (a.active()) Defect<tags::GPU::CUDA>::value(r.vector(), rhs.vector(), a.inner_matrix(), b.vector());

    MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);
    recv_requests.clear();
    DenseVector<DT_> missing_values(a.outer_matrix().columns());
    void * missing_values_gpu(missing_values.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    cuda_upload(missing_values_array, missing_values_gpu, missing_values.size() * sizeof(DT_));
    missing_values.unlock(lm_write_only);

    // berechne aeussere anteile
    DenseVector<DT_> r_outer(r.local_size());
    if (a.active()) Product<tags::GPU::CUDA>::value(r_outer, a.outer_matrix(), missing_values);
    if (a.active()) Difference<tags::GPU::CUDA>::value(r.vector(), r_outer);

    MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
    send_requests.clear();
    cuda_free_host(missing_values_array);
    cuda_free_host(b_cpu);
    delete[] send_data;
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


template struct MPIOps<tags::CPU>;
template void MPIOps<tags::CPU>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

template struct MPIOps<tags::CPU::SSE>;
template void MPIOps<tags::CPU::SSE>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::SSE>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::SSE>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::SSE>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::SSE>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::SSE>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::SSE>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::SSE>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

template struct MPIOps<tags::CPU::MultiCore::SSE>;
template void MPIOps<tags::CPU::MultiCore::SSE>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::SSE>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::CPU::MultiCore::SSE>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::MultiCore::SSE>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::MultiCore::SSE>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::CPU::MultiCore::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::MultiCore::SSE>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::MultiCore::SSE>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

#ifdef HONEI_CUDA
template struct MPIOps<tags::GPU::CUDA>;
template void MPIOps<tags::GPU::CUDA>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::GPU::CUDA>::defect(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & rhs, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template double MPIOps<tags::GPU::CUDA>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::GPU::CUDA>::element_product(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::GPU::CUDA>::norm_l2_false(const DenseVectorMPI<double> & x);
template void MPIOps<tags::GPU::CUDA>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::GPU::CUDA>::scale(DenseVectorMPI<double> & x, double a);
template void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::GPU::CUDA>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::GPU::CUDA>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
#endif
