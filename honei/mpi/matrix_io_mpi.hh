/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef MPI_GUARD_MATRIX_IO_MPI_HH
#define MPI_GUARD_MATRIX_IO_MPI_HH 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <honei/math/matrix_io.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>
#include <honei/mpi/sparse_matrix_csr_mpi.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/algorithm.hh>
#include <vector>
#include <honei/backends/mpi/operations.hh>

using namespace honei;

template<typename IOFormat_>
class MatrixIOMPI
{
};

template<>
class MatrixIOMPI<io_formats::ELL>
{
    public:
        static unsigned long read_matrix_rows(std::string input)
        {
            unsigned long _rank(mpi::mpi_comm_rank());
            unsigned long rows(0);
            if (_rank == 0)
            {
                FILE* file(NULL);
                file = fopen(input.c_str(), "rb");
                if (file == NULL)
                    throw InternalError("File "+input+" not found!");
                uint64_t size;
                uint64_t rowst;
                int status = fread(&size, sizeof(uint64_t), 1, file);
                status = fread(&rowst, sizeof(uint64_t), 1, file);
                fclose(file);
                rows = rowst;
                mpi::mpi_bcast(&rows, 1, 0);
            }
            else
            {
                mpi::mpi_bcast(&rows, 1, 0);
            }
            return rows;
        }

        template <typename DT_>
            static SparseMatrix<DT_> read_matrix(std::string input, HONEI_UNUSED DT_ datatype)
            {
                unsigned long _rank(mpi::mpi_comm_rank());
                if (_rank == 0)
                {
                    FILE* file(NULL);
                    file = fopen(input.c_str(), "rb");
                    if (file == NULL)
                        throw InternalError("File "+input+" not found!");
                    uint64_t size;
                    uint64_t rows;
                    uint64_t columns;
                    uint64_t stride;
                    uint64_t num_cols_per_row;
                    int status = fread(&size, sizeof(uint64_t), 1, file);
                    status = fread(&rows, sizeof(uint64_t), 1, file);
                    status = fread(&columns, sizeof(uint64_t), 1, file);
                    status = fread(&stride, sizeof(uint64_t), 1, file);
                    status = fread(&num_cols_per_row, sizeof(uint64_t), 1, file);
                    long int pos_aj(ftell(file));
                    fseek(file, size * sizeof(uint64_t), SEEK_CUR);
                    long int pos_ax(ftell(file));

                    unsigned long local_rows(DenseVectorMPI<DT_>::calc_size(rows));
                    unsigned long local_offset(DenseVectorMPI<DT_>::calc_offset(rows));
                    unsigned long alien_rows(rows);
                    unsigned long alien_columns(columns);
                    unsigned long alien_stride(stride);
                    unsigned long alien_num_cols_per_row(num_cols_per_row);

                    mpi::mpi_bcast(&alien_rows, 1, 0);
                    mpi::mpi_bcast(&alien_columns, 1, 0);
                    mpi::mpi_bcast(&alien_stride, 1, 0);
                    mpi::mpi_bcast(&alien_num_cols_per_row, 1, 0);

                    SparseMatrix<DT_> local_matrix(local_rows, columns);

                    uint64_t * indices = new uint64_t[stride];
                    double * values = new double[stride];
                    for (unsigned long current_col(0) ; current_col < num_cols_per_row ; ++current_col)
                    {
                        fseek(file, pos_aj, SEEK_SET);
                        fseek(file, (current_col * stride) * sizeof(uint64_t), SEEK_CUR);

                        int status = fread(indices, sizeof(uint64_t), stride, file);

                        fseek(file, pos_ax, SEEK_SET);
                        fseek(file, (current_col * stride) * sizeof(double), SEEK_CUR);

                        status = fread(values, sizeof(double), stride, file);

                        mpi::mpi_bcast(indices, stride, _rank);
                        mpi::mpi_bcast(values, stride, _rank);

                        for (unsigned long i(local_offset) ; i < local_offset + local_rows ; ++i)
                        {
                            if (values[i] != double(0))
                                local_matrix(i - local_offset, indices[i], values[i]);
                        }
                    }
                    delete[] indices;
                    delete[] values;

                    fclose(file);
                    return local_matrix;
                }
                else
                {
                    unsigned long rows;
                    unsigned long columns;
                    unsigned long stride;
                    unsigned long num_cols_per_row;
                    mpi::mpi_bcast(&rows, 1, 0);
                    mpi::mpi_bcast(&columns, 1, 0);
                    mpi::mpi_bcast(&stride, 1, 0);
                    mpi::mpi_bcast(&num_cols_per_row, 1, 0);
                    unsigned long local_rows(DenseVectorMPI<DT_>::calc_size(rows));
                    unsigned long local_offset(DenseVectorMPI<DT_>::calc_offset(rows));
                    SparseMatrix<DT_> local_matrix(local_rows, columns);

                    uint64_t * indices = new uint64_t[stride];
                    double * values = new double[stride];

                    for (unsigned long current_col(0) ; current_col < num_cols_per_row ; ++current_col)
                    {
                        mpi::mpi_bcast(indices, stride, 0);
                        mpi::mpi_bcast(values, stride, 0);

                        for (unsigned long i(local_offset) ; i < local_offset + local_rows ; ++i)
                        {
                            if (values[i] != double(0))
                                local_matrix(i - local_offset, indices[i], values[i]);
                        }
                    }
                    delete[] indices;
                    delete[] values;

                    return local_matrix;
                }
            }
};

template<typename IOFormat_>
class MatrixIOMPI_ELL
{
};

template<>
class MatrixIOMPI_ELL<io_formats::ELL>
{
    public:
        template <typename DT_>
            static SparseMatrixELLMPI<DT_> read_matrix(std::string input, HONEI_UNUSED DT_ datatype)
            {
                unsigned long global_rows = MatrixIOMPI<io_formats::ELL>::read_matrix_rows(input);
                SparseMatrix<DT_> local = MatrixIOMPI<io_formats::ELL>::read_matrix(input, datatype);
                SparseMatrixELLMPI<DT_> result(local, global_rows);
                return result;
            }
};

template<typename IOFormat_>
class MatrixIOMPI_CSR
{
};

template<>
class MatrixIOMPI_CSR<io_formats::ELL>
{
    public:
        template <typename DT_>
            static SparseMatrixCSRMPI<DT_> read_matrix(std::string input, HONEI_UNUSED DT_ datatype)
            {
                unsigned long global_rows = MatrixIOMPI<io_formats::ELL>::read_matrix_rows(input);
                SparseMatrix<DT_> local = MatrixIOMPI<io_formats::ELL>::read_matrix(input, datatype);
                SparseMatrixCSRMPI<DT_> result(local, global_rows);
                return result;
            }
};
#endif
