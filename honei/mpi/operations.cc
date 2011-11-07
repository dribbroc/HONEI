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
#include <honei/backends/mpi/operations.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/product.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>

using namespace honei;

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::difference(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Difference<Tag_>::value(r.vector(), x.vector(), y.vector());
}

template <typename Tag_>
template <typename DT_>
DT_ MPIOps<Tag_>::dot_product(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    DT_ local_result(DotProduct<Tag_>::value(x.vector(), y.vector()));
    DT_ result(DT_(0));

    MPI_Allreduce(&local_result, &result, 1, mpi::MPIType<DT_>::value(), MPI_SUM, MPI_COMM_WORLD);

    return result;
}

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::product(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
{
    // cast away const as mpi does not like const pointers
    DT_ * bp = (DT_*)b.elements();

    int myrank;
    mpi::mpi_comm_rank(&myrank);
    int com_size;
    mpi::mpi_comm_size(&com_size);

    std::vector<MPI_Request> requests;

    // \TODO use only one send/recv per process
    DT_ missing_values[a.missing_indices().size()];
    // empfange alle fehlenden werte
    {
        unsigned long j(0);
        for (std::set<unsigned long>::iterator i(a.missing_indices().begin()) ; i != a.missing_indices().end() ; ++i, ++j)
        {
            requests.push_back(mpi::mpi_irecv(missing_values + j, 1, MPI_ANY_SOURCE, *i));
        }
    }

    // sende alle werte, die anderen fehlen
    for (int rank(0) ; rank < com_size ; ++rank)
    {
        if (rank == myrank)
            continue;

        unsigned long j(0);
        for (std::set<unsigned long>::iterator i(a.alien_indices(rank).begin()) ; i != a.alien_indices(rank).end() ; ++i, ++j)
        {
            requests.push_back(mpi::mpi_isend(bp + *i - a.offset(), 1, rank, *i));
        }
    }

    // berechne innere anteile
    Product<Tag_>::value(r.vector(), a.inner_matrix(), b.vector(), true);

    // TODO nur auf empfang warten - senden warten reicht auch wenn ich ganz fertig bin.
    MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

    // berechne aeussere anteile
    unsigned long ix(0);
    DT_ * r_ele(r.elements());
    const DT_ * b_ele(b.elements());
    //vor eigenem vektor
    for (unsigned long coli(0) ; coli < a.before_cols().size() ; ++coli)
    {
        unsigned long col(a.before_cols().at(coli));
        const SparseVector<DT_> scol((a.outer_matrix()).column(col));
        const unsigned long ue(scol.used_elements());
        const DT_* scol_ele(scol.elements());
        const unsigned long * indices(scol.indices());

        for (unsigned long row(0) ; row < ue ; ++row)
        {
            const unsigned long row_index(indices[row]);
            r_ele[row_index] +=  scol_ele[row] * missing_values[ix];
        }
        ++ix;
    }

    //in eigenem vektor
    for (unsigned long coli(0) ; coli < a.middle_cols().size() ; ++coli)
    {
        unsigned long col(a.middle_cols().at(coli));
        const SparseVector<DT_> scol((a.outer_matrix()).column(col));
        const unsigned long ue(scol.used_elements());
        const unsigned long * indices(scol.indices());
        const DT_* scol_ele(scol.elements());

        for (unsigned long row(0) ; row < ue ; ++row)
        {
            const unsigned long row_index(indices[row]);
            r_ele[row_index] += scol_ele[row] * b_ele[ col  - a.offset()];
        }
    }

    //nach eigenem vektor
    for (unsigned long coli(0) ; coli < a.behind_cols().size() ; ++coli)
    {
        unsigned long col(a.behind_cols().at(coli));
        const SparseVector<DT_> scol((a.outer_matrix()).column(col));
        const unsigned long ue(scol.used_elements());
        const unsigned long * indices(scol.indices());
        const DT_* scol_ele(scol.elements());

        for (unsigned long row(0) ; row < ue ; ++row)
        {
            const unsigned long row_index(indices[row]);
            r_ele[row_index] +=  scol_ele[row] * missing_values[ix];
        }
        ++ix;
    }
}

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::scaled_sum(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a)
{
    ScaledSum<Tag_>::value(r.vector(), x.vector(), y.vector(), a);
}

template <typename Tag_>
template <typename DT_>
void MPIOps<Tag_>::sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
{
    Sum<Tag_>::value(x.vector(), y.vector());
}

// \ TODO use atomic defect operation

template struct MPIOps<tags::CPU>;
template void MPIOps<tags::CPU>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);


template struct MPIOps<tags::CPU::SSE>;
template void MPIOps<tags::CPU::SSE>::difference(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template double MPIOps<tags::CPU::SSE>::dot_product(const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);
template void MPIOps<tags::CPU::SSE>::product(DenseVectorMPI<double> & r, const SparseMatrixELLMPI<double> & a, const DenseVectorMPI<double> & b);
template void MPIOps<tags::CPU::SSE>::scaled_sum(DenseVectorMPI<double> & r, const DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y, double a);
template void MPIOps<tags::CPU::SSE>::sum(DenseVectorMPI<double> & x, const DenseVectorMPI<double> & y);

