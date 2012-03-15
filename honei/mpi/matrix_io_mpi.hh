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
            FILE* file(NULL);
            file = fopen(input.c_str(), "rb");
            if (file == NULL)
                throw InternalError("File "+input+" not found!");
            uint64_t size;
            uint64_t rows;
            int status = fread(&size, sizeof(uint64_t), 1, file);
            status = fread(&rows, sizeof(uint64_t), 1, file);
            fclose(file);
            return rows;
        }

        template <typename DT_>
            static SparseMatrix<DT_> read_matrix(std::string input, HONEI_UNUSED DT_ datatype)
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

                SparseMatrix<DT_> local_matrix(local_rows, columns);
                for (unsigned long current_row(local_offset) ; current_row < local_rows + local_offset ; ++current_row)
                {
                    for (unsigned long cols(0) ; cols < num_cols_per_row ; ++cols)
                    {
                        fseek(file, pos_aj, SEEK_SET);
                        fseek(file, (current_row + (cols * stride)) * sizeof(uint64_t), SEEK_CUR);

                        uint64_t temp;
                        int status = fread(&temp, sizeof(uint64_t), 1, file);
                        unsigned long icol(temp);

                        fseek(file, pos_ax, SEEK_SET);
                        fseek(file, (current_row + (cols * stride)) * sizeof(double), SEEK_CUR);

                        double ival;
                        status = fread(&ival, sizeof(double), 1, file);
                        if (ival != DT_(0)) local_matrix(current_row-local_offset, icol, ival);
                    }
                }

                fclose(file);
                return local_matrix;
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
                unsigned long global_rows(MatrixIOMPI<io_formats::ELL>::read_matrix_rows(input));
                SparseMatrix<DT_> local(MatrixIOMPI<io_formats::ELL>::read_matrix(input, datatype));
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
                unsigned long global_rows(MatrixIOMPI<io_formats::ELL>::read_matrix_rows(input));
                SparseMatrix<DT_> local(MatrixIOMPI<io_formats::ELL>::read_matrix(input, datatype));
                SparseMatrixCSRMPI<DT_> result(local, global_rows);
                return result;
            }
};
#endif
