/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#pragma once
#ifndef MPI_GUARD_OPERATIONS_HH
#define MPI_GUARD_OPERATIONS_HH 1

#include <mpi.h>

namespace honei
{
    namespace mpi
    {
        void mpi_init(int * argc, char*** argv);
        void mpi_init();
        void mpi_finalize();
        void mpi_comm_size(int * size, MPI_Comm com = MPI_COMM_WORLD);
        void mpi_comm_rank(int * id, MPI_Comm com = MPI_COMM_WORLD);
        template <typename DT_> void mpi_bcast(DT_ * data, unsigned long size, int sender, MPI_Comm com = MPI_COMM_WORLD);
        template <typename DT_> void mpi_send(DT_ * data, unsigned long size, int target, int tag, MPI_Comm com = MPI_COMM_WORLD);
        template <typename DT_> void mpi_recv(DT_ * data, unsigned long size, int sender, int tag, MPI_Comm com = MPI_COMM_WORLD);
        template <typename DT_> MPI_Request mpi_isend(DT_ * data, unsigned long size, int target, int tag, MPI_Comm com = MPI_COMM_WORLD);
        template <typename DT_> MPI_Request mpi_irecv(DT_ * data, unsigned long size, int sender, int tag, MPI_Comm com = MPI_COMM_WORLD);

        template <typename DT_>
        class MPIType
        {
        };

        template <>
        class MPIType<float>
        {
            public:
                static inline MPI_Datatype value()
                {
                    return MPI_FLOAT;
                }
        };

        template <>
        class MPIType<double>
        {
            public:
                static inline MPI_Datatype value()
                {
                    return MPI_DOUBLE;
                }
        };

        template <>
        class MPIType<unsigned long>
        {
            public:
                static inline MPI_Datatype value()
                {
                    return MPI_UNSIGNED_LONG;
                }
        };

    }
}

#endif
