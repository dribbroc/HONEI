/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <mpi.h>
#include <honei/backends/mpi/operations.hh>
#ifdef HONEI_GMP
#include <gmpxx.h>
#endif

namespace honei
{
    namespace mpi
    {
        void mpi_init(int * argc, char*** argv)
        {
            MPI_Init(argc, argv);
        }

        void mpi_init()
        {
            int argc(0);
            char argv = 'a';
            char * argv2(&argv);
            char ** argv3(&argv2);

            MPI_Init(&argc, &argv3);
        }

        void mpi_finalize()
        {
            MPI_Finalize();
        }

        void mpi_comm_size(int * size, MPI_Comm com)
        {
            MPI_Comm_size(com, size);
        }

        int mpi_comm_size(MPI_Comm com)
        {
            int size;
            MPI_Comm_size(com, &size);
            return size;
        }

        void mpi_comm_rank(int * id, MPI_Comm com)
        {
            MPI_Comm_rank(com, id);
        }

        int mpi_comm_rank(MPI_Comm com)
        {
            int id;
            MPI_Comm_rank(com, &id);
            return id;
        }

        template <typename DT_>
        void mpi_bcast(DT_ * data, unsigned long size, int sender, MPI_Comm com)
        {
            MPI_Bcast(data, size, MPIType<DT_>::value(), sender, com);
        }

#ifdef HONEI_GMP
        template <>
        void mpi_bcast(mpf_class * data, unsigned long size, int sender, MPI_Comm com)
        {
            double tdata[size];
            if (mpi_comm_rank() == sender)
                for (unsigned long i(0) ; i < size ; ++i)
                    tdata[i] = mpf_get_d(data[i].get_mpf_t());

            MPI_Bcast(tdata, size, MPIType<double>::value(), sender, com);

            if (mpi_comm_rank() != sender)
                for (unsigned long i(0) ; i < size ; ++i)
                    data[i] = tdata[i];
        }
#endif

        template <typename DT_>
        void mpi_send(DT_ * data, unsigned long size, int target, int tag, MPI_Comm com)
        {
            MPI_Send(data, size, MPIType<DT_>::value(), target, tag, com);
        }

        template <typename DT_>
        void mpi_recv(DT_ * data, unsigned long size, int sender, int tag, MPI_Comm com)
        {
            MPI_Status _stat;
            MPI_Recv(data, size, MPIType<DT_>::value(), sender, tag, com, &_stat);
        }

        template <typename DT_> MPI_Request mpi_isend(DT_ * data, unsigned long size, int target, int tag, MPI_Comm com)
        {
            MPI_Request request;
            MPI_Isend(data, size, mpi::MPIType<DT_>::value(), target, tag, com, &(request));
            return request;
        }

        template <typename DT_> MPI_Request mpi_irecv(DT_ * data, unsigned long size, int sender, int tag, MPI_Comm com)
        {
            MPI_Request request;
            MPI_Irecv(data, size, mpi::MPIType<DT_>::value(), sender, tag, com, &(request));
            return request;
        }

        template void mpi_bcast<float>(float * data, unsigned long size, int sender, MPI_Comm com);
        template void mpi_bcast<double>(double * data, unsigned long size, int sender, MPI_Comm com);
        template void mpi_bcast<unsigned long>(unsigned long * data, unsigned long size, int sender, MPI_Comm com);

        template void mpi_send<float>(float * data, unsigned long size, int target, int tag, MPI_Comm com);
        template void mpi_send<double>(double * data, unsigned long size, int target, int tag, MPI_Comm com);
        template void mpi_send<unsigned long>(unsigned long * data, unsigned long size, int target, int tag, MPI_Comm com);

        template void mpi_recv<float>(float * data, unsigned long size, int sender, int tag, MPI_Comm com);
        template void mpi_recv<double>(double * data, unsigned long size, int sender, int tag, MPI_Comm com);
        template void mpi_recv<unsigned long>(unsigned long * data, unsigned long size, int sender, int tag, MPI_Comm com);

        template MPI_Request mpi_isend<float>(float * data, unsigned long size, int target, int tag, MPI_Comm com);
        template MPI_Request mpi_isend<double>(double * data, unsigned long size, int target, int tag, MPI_Comm com);
        template MPI_Request mpi_isend<unsigned long>(unsigned long * data, unsigned long size, int target, int tag, MPI_Comm com);

        template MPI_Request mpi_irecv<float>(float * data, unsigned long size, int sender, int tag, MPI_Comm com);
        template MPI_Request mpi_irecv<double>(double * data, unsigned long size, int sender, int tag, MPI_Comm com);
        template MPI_Request mpi_irecv<unsigned long>(unsigned long * data, unsigned long size, int sender, int tag, MPI_Comm com);
    }
}
