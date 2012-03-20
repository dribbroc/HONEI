/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef MATH_GUARD_MATRIX_IO_HH
#define MATH_GUARD_MATRIX_IO_HH 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/algorithm.hh>
#include <honei/util/attributes.hh>
#include <honei/util/configuration.hh>
#include <vector>
#include <algorithm>
#ifdef HONEI_MPI
#include <honei/backends/mpi/operations.hh>
#endif

using namespace honei;

namespace io_formats
{
    class MTX;
    class M;
    class ELL;
}

template<typename IOFormat_>
class MatrixIO
{
};

template<>
class MatrixIO<io_formats::M>
{
    public:
        template<typename DT_>
            static SparseMatrix<DT_> read_matrix(std::string filename, DT_)
            {
                std::ifstream file(filename.c_str());
                if (!file.is_open())
                    throw honei::InternalError("Unable to open MATLAB file: " + filename);

                std::vector<unsigned long> row_indices;
                std::vector<unsigned long> column_indices;
                std::vector<DT_> data;

                ///Attention: MatrixMarket indices are 1-based!!!
                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);

                    if(line.find("]", 0) < line.npos)
                        break;

                    if(line.find("data", 0) < line.npos)
                    {
                        continue;
                    }

                    if(file.eof())
                        break;

                    std::string c_s, r_s, n_z_s;

                    std::string::size_type first_digit(line.find_first_not_of(" "));
                    line.erase(0, first_digit);
                    std::string::size_type first_blank(line.find_first_of(" "));
                    for(unsigned long i(0) ; i < first_blank ; ++i)
                    {
                        r_s.append(1, line[i]);
                    }
                    line.erase(0, first_blank + 1);
                    first_digit = line.find_first_not_of(" ");
                    line.erase(0, first_digit);

                    std::string::size_type second_blank(line.find_first_of(" "));
                    for(unsigned long i(0) ; i < second_blank ; ++i)
                    {
                        c_s.append(1, line[i]);
                    }
                    line.erase(0, second_blank + 1);
                    first_digit = line.find_first_not_of(" ");
                    line.erase(0, first_digit);


                    std::string::size_type first_semicolon(line.find_first_of(";"));
                    for(unsigned long i(0) ; i < first_semicolon ; ++i)
                    {
                        n_z_s.append(1, line[i]);
                    }


                    unsigned long c = (unsigned long)atol(c_s.c_str());
                    unsigned long r = (unsigned long)atol(r_s.c_str());
                    DT_ n_z = (DT_)atof(n_z_s.c_str());

                    row_indices.push_back(r-1);
                    column_indices.push_back(c-1);
                    data.push_back(n_z);

                }
                file.close();
                unsigned long columns(*std::max_element(column_indices.begin(), column_indices.end()) + 1);
                unsigned long rows(*std::max_element(row_indices.begin(), row_indices.end()) + 1);
                SparseMatrix<DT_> result(rows, columns, &row_indices[0], &column_indices[0], &data[0], data.size());
                return result;
            }

};

//MATRIX MARKET TYPE
template<>
class MatrixIO<io_formats::MTX>
{
    private:
        static void get_sizes(std::string filename, unsigned long & r,
                unsigned long & c,
                unsigned long & n_z,
                unsigned long & data_begin)
        {
            std::string c_s, r_s, n_z_s;
            std::ifstream file(filename.c_str());

            unsigned long data_index(0);

            if(file.is_open())
            {
                while(true)
                {
                    std::string line;
                    std::getline(file, line);
                    if(line.find("%", 0) < line.npos)
                    {
                        ++data_index;
                        continue;
                    }
                    else
                    {
                        std::string::size_type first_blank(line.find_first_of(" "));
                        for(unsigned long i(0) ; i < first_blank ; ++i)
                        {
                            r_s.append(1, line[i]);
                        }
                        line.erase(0, first_blank + 1);

                        std::string::size_type second_blank(line.find_first_of(" "));
                        for(unsigned long i(0) ; i < second_blank ; ++i)
                        {
                            c_s.append(1, line[i]);
                        }
                        line.erase(0, second_blank + 1);

                        std::string::size_type eol(line.length());
                        for(unsigned long i(0) ; i < eol ; ++i)
                        {
                            n_z_s.append(1, line[i]);
                        }

                        ++data_index;
                        break;
                    }
                }

                file.close();

                c = (unsigned long)atol(c_s.c_str());
                r = (unsigned long)atol(r_s.c_str());
                n_z = (unsigned long)atol(n_z_s.c_str());

                data_begin = data_index;
            }
            else
                throw honei::InternalError("Unable to open MatrixMarket file.");
        }

    public:
        ///Read in sparse data only
        template<typename DT_>
            static SparseMatrix<DT_> read_matrix(std::string filename, DT_)
            {
                unsigned long cols(0), rows(0);

                std::ifstream file(filename.c_str());
                if (!file.is_open())
                    throw honei::InternalError("Unable to open Matrixmarket M file: "+ filename);

                ///Attention: MatrixMarket indices are 1-based!!!

                std::vector<unsigned long> column_indices;
                std::vector<unsigned long> row_indices;
                std::vector<DT_> data;
                bool first(true);
                bool symmetric(false);

                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);
                    if(line.find("%", 0) < line.npos)
                    {
                        if (line.find("symmetric") < line.npos)
                        {
                            symmetric = true;
                        }
                        continue;
                    }
                    if(file.eof())
                        break;

                    if (first)
                    {
                        std::string c_s, r_s;
                        std::string::size_type first_blank(line.find_first_of(" "));
                        for(unsigned long i(0) ; i < first_blank ; ++i)
                        {
                            r_s.append(1, line[i]);
                        }
                        line.erase(0, first_blank + 1);

                        std::string::size_type second_blank(line.find_first_of(" "));
                        for(unsigned long i(0) ; i < second_blank ; ++i)
                        {
                            c_s.append(1, line[i]);
                        }

                        cols = (unsigned long)atol(c_s.c_str());
                        rows = (unsigned long)atol(r_s.c_str());
                        first = false;
                        continue;
                    }

                    else
                    {
                        std::string c_s, r_s, n_z_s;

                        std::string::size_type first_blank(line.find_first_of(" "));
                        for(unsigned long j(0) ; j < first_blank ; ++j)
                        {
                            r_s.append(1, line[j]);
                        }
                        line.erase(0, first_blank + 1);

                        std::string::size_type second_blank(line.find_first_of(" "));
                        for(unsigned long j(0) ; j < second_blank ; ++j)
                        {
                            c_s.append(1, line[j]);
                        }
                        line.erase(0, second_blank + 1);

                        std::string::size_type eol(line.length());
                        for(unsigned long j(0) ; j < eol ; ++j)
                        {
                            n_z_s.append(1, line[j]);
                        }


                        unsigned long c = (unsigned long)atol(c_s.c_str());
                        unsigned long r = (unsigned long)atol(r_s.c_str());
                        DT_ n_z = (DT_)atof(n_z_s.c_str());

                        row_indices.push_back(r - 1);
                        column_indices.push_back(c -1);
                        data.push_back(n_z);
                    }

                }
                file.close();

                if (symmetric)
                {
                    const unsigned long old_size(data.size());
                    for (unsigned long i(0) ; i < old_size ; ++i)
                    {
                        if (! (row_indices.at(i) == column_indices.at(i)))
                        {
                            row_indices.push_back(column_indices.at(i));
                            column_indices.push_back(row_indices.at(i));
                            data.push_back(data.at(i));
                        }
                    }
                }

                SparseMatrix<DT_> result(rows, cols, &row_indices[0], &column_indices[0], &data[0], data.size());
                return result;
            }


        template<typename DT_>
            static DenseMatrix<DT_> read_matrix(std::string filename, DT_ base, unsigned long & non_zeros)
            {
                unsigned long columns, rows, non_data_lines;
                get_sizes(filename, rows, columns, non_zeros, non_data_lines);

                DenseMatrix<DT_> result(rows, columns, base);

                std::ifstream file(filename.c_str());
                file.is_open();

                for(unsigned long i(0) ; i < non_data_lines ; ++i)
                {
                    std::string line;
                    std::getline(file, line);
                }

                ///Attention: MatrixMarket indices are 1-based!!!
                for(unsigned long i(0) ; i < non_zeros ; ++i)
                {
                    std::string line;
                    std::getline(file, line);

                    std::string c_s, r_s, n_z_s;

                    std::string::size_type first_blank(line.find_first_of(" "));
                    for(unsigned long j(0) ; j < first_blank ; ++j)
                    {
                        r_s.append(1, line[j]);
                    }
                    line.erase(0, first_blank + 1);

                    std::string::size_type second_blank(line.find_first_of(" "));
                    for(unsigned long j(0) ; j < second_blank ; ++j)
                    {
                        c_s.append(1, line[j]);
                    }
                    line.erase(0, second_blank + 1);

                    std::string::size_type eol(line.length());
                    for(unsigned long j(0) ; j < eol ; ++j)
                    {
                        n_z_s.append(1, line[j]);
                    }


                    unsigned long c = (unsigned long)atol(c_s.c_str());
                    unsigned long r = (unsigned long)atol(r_s.c_str());
                    DT_ n_z = (DT_)atof(n_z_s.c_str());

                    result[r - 1][c - 1] = n_z;
                }
                file.close();

                return result;
            }

        template<typename DT_>
            static void write_matrix(std::string filename, const SparseMatrix<DT_> & matrix)
            {
                unsigned long used_elements(0);
                for (unsigned long i(0) ; i < matrix.rows() ; ++i)
                    used_elements += matrix[i].used_elements();

                FILE* file;
                file = fopen(filename.c_str(), "w");
                std::string header("%%MatrixMarket matrix coordinate real general\n");
                std::string mnl;
                mnl += stringify(matrix.rows());
                mnl += " ";
                mnl += stringify(matrix.columns());
                mnl += " ";
                mnl += stringify(used_elements);
                mnl += "\n";
                fprintf(file, "%s", header.c_str());
                fprintf(file, "%s", mnl.c_str());
                for (typename SparseMatrix<DT_>::NonZeroConstElementIterator i(matrix.begin_non_zero_elements()), i_end(matrix.end_non_zero_elements()) ;
                        i != i_end ; ++i)
                {
                    std::string temp;
                    temp += stringify(i.row() + 1) + " ";
                    temp += stringify(i.column() + 1) + " ";
                    temp += stringify(*i) + "\n";
                    fprintf(file, "%s", temp.c_str());
                }
                fclose(file);
            }
};

template<>
class MatrixIO<io_formats::ELL>
{
    public:
        template <typename DT_>
            static void write_matrix(std::string & output, SparseMatrixELL<DT_> & smatrix)
            {
                if (sizeof(DT_) != 8)
                    throw InternalError("Only double ell output supported!");
                else if (sizeof(unsigned long) != 8)
                    throw InternalError("Only 64 bit machine output supported!");
                /// \todo convert uint32 to uint64 if needed (see read_matrix for conversion)
                else if (smatrix.threads() != 1)
                    throw InternalError("Only Matrices with 1 threads data layout are supported for export");

                else
                {
                    FILE* file;
                    file = fopen(output.c_str(), "wb");
                    uint64_t size(smatrix.Aj().size());
                    uint64_t rows(smatrix.rows());
                    uint64_t columns(smatrix.columns());
                    uint64_t stride(smatrix.stride());
                    uint64_t num_cols_per_row(smatrix.num_cols_per_row());
                    fwrite(&size, sizeof(uint64_t), 1, file);
                    fwrite(&rows, sizeof(uint64_t), 1, file);
                    fwrite(&columns, sizeof(uint64_t), 1, file);
                    fwrite(&stride, sizeof(uint64_t), 1, file);
                    fwrite(&num_cols_per_row, sizeof(uint64_t), 1, file);
                    fwrite(smatrix.Aj().elements(), sizeof(uint64_t), size, file);
                    fwrite(smatrix.Ax().elements(), sizeof(double), size, file);
                    fclose(file);
                }
            }

        template <typename DT_>
            static SparseMatrixELL<DT_> read_matrix(std::string input, HONEI_UNUSED DT_ datatype)
            {
#ifdef HONEI_MPI
                int rank(mpi::mpi_comm_rank());
                if(rank == 0)
                {
#endif
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
                    DenseVector<unsigned long> ajc(size);
                    if (sizeof(unsigned long) == sizeof(uint64_t))
                    {
                        status = fread(ajc.elements(), sizeof(uint64_t), size, file);
                    }
                    else
                    {
                        uint64_t aj[size];
                        status = fread(aj, sizeof(uint64_t), size, file);
                        for (unsigned long i(0) ; i < size ; ++i)
                        {
                            ajc[i] = aj[i];
                        }
                    }
                    DenseVector<double> ax(size);
                    status = fread(ax.elements(), sizeof(double), size, file);
                    fclose(file);
                    DenseVector<DT_> axc(size);
                    unsigned long crows(rows);
                    unsigned long ccolumns(columns);
                    unsigned long cstride(stride);
                    unsigned long cnum_cols_per_row(num_cols_per_row);
                    convert<tags::CPU>(axc, ax);
#ifdef HONEI_MPI
                    unsigned long vec_size(ajc.size());
                    mpi::mpi_bcast(&vec_size, 1, 0);
                    mpi::mpi_bcast(&crows, 1, 0);
                    mpi::mpi_bcast(&ccolumns, 1, 0);
                    mpi::mpi_bcast(&cstride, 1, 0);
                    mpi::mpi_bcast(&cnum_cols_per_row, 1, 0);

                    mpi::mpi_bcast(ajc.elements(), vec_size, 0);
                    mpi::mpi_bcast(axc.elements(), vec_size, 0);

                    SparseMatrixELL<DT_> smatrix(crows, ccolumns, cstride, cnum_cols_per_row, ajc, axc, 1);
                    if (Configuration::instance()->get_value("ell::threads", 1) != 1)
                    {
                        SparseMatrix<DT_> bla(smatrix);
                        SparseMatrixELL<DT_> smatrix2(bla);
                        return smatrix2;
                    }
                    return smatrix;
                }
                else
                {
                    unsigned long vec_size, crows, ccolumns, cstride, cnum_cols_per_row;
                    mpi::mpi_bcast(&vec_size, 1, 0);
                    mpi::mpi_bcast(&crows, 1, 0);
                    mpi::mpi_bcast(&ccolumns, 1, 0);
                    mpi::mpi_bcast(&cstride, 1, 0);
                    mpi::mpi_bcast(&cnum_cols_per_row, 1, 0);

                    DenseVector<unsigned long> ajc(vec_size);
                    DenseVector<DT_> axc(vec_size);

                    mpi::mpi_bcast(ajc.elements(), vec_size, 0);
                    mpi::mpi_bcast(axc.elements(), vec_size, 0);

                    SparseMatrixELL<DT_> smatrix(crows, ccolumns, cstride, cnum_cols_per_row, ajc, axc, 1);
                    if (Configuration::instance()->get_value("ell::threads", 1) != 1)
                    {
                        SparseMatrix<DT_> bla(smatrix);
                        SparseMatrixELL<DT_> smatrix2(bla);
                        return smatrix2;
                    }
                    return smatrix;
                }
#endif
#ifndef HONEI_MPI
                SparseMatrixELL<DT_> smatrix(crows, ccolumns, cstride, cnum_cols_per_row, ajc, axc, 1);
                if (Configuration::instance()->get_value("ell::threads", 1) != 1)
                {
                    SparseMatrix<DT_> bla(smatrix);
                    SparseMatrixELL<DT_> smatrix2(bla);
                    return smatrix2;
                }
                return smatrix;
#endif
            }
};
#endif
