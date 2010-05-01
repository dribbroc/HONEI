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
        static void get_sizes(std::string filename, unsigned long & r,
                                                    unsigned long & c,
                                                    unsigned long & n_z,
                                                    unsigned long & data_begin)
        {
            n_z = 0;
            data_begin = 0;

            std::string c_s, r_s;
            std::ifstream file(filename.c_str());


            if(file.is_open())
            {
                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);

                    //ignore header and footer
                    if(line.find("]", 0) < line.npos)
                        break;

                    if(line.find("data", 0) < line.npos)
                    {
                        ++data_begin;
                        continue;
                    }
                    else
                    {
                        ++n_z;
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

                        c = (unsigned long)atol(c_s.c_str());
                        r = (unsigned long)atol(r_s.c_str());

                        c_s.clear();
                        r_s.clear();

                    }

                }

                file.close();

            }
            else
                throw honei::InternalError("Unable to open MATLAB file.");
        }

        static unsigned long get_non_zeros(std::string filename)
        {
            unsigned long n_z(0);
            std::ifstream file(filename.c_str());
            if(file.is_open())
            {
                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);

                    //ignore header and footer
                    if(line.find("]", 0) < line.npos)
                        break;

                    if(line.find("data", 0) < line.npos)
                    {
                        continue;
                    }
                    else
                    {
                        ++n_z;
                    }

                }

                file.close();

            }
            else
                throw honei::InternalError("Unable to open MATLAB file.");

            return n_z;
        }

        template<typename DT_>
            static void read_matrix(std::string filename, DenseVector<unsigned long> & row_indices,
                                                          DenseVector<unsigned long> & column_indices,
                                                          DenseVector<DT_> & data)
            {
                unsigned long columns, rows, non_data_lines, non_zeros;
                get_sizes(filename, rows, columns, non_zeros, non_data_lines);

                std::ifstream file(filename.c_str());
                file.is_open();

                for(unsigned long i(0) ; i < non_data_lines ; ++i)
                {
                    std::string line;
                    std::getline(file, line);
                }

                ///Attention: MatrixMarket indices are 1-based!!!
                for(unsigned long j(0) ; j < non_zeros ; ++j)
                {
                    std::string line;
                    std::getline(file, line);

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

                    row_indices[j] = r - 1;
                    column_indices[j] = c -1;
                    data[j] = n_z;

                }
                file.close();
            }
};

//MATRIX MARKET TYPE
template<>
class MatrixIO<io_formats::MTX>
{
    public:
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

        ///Only read size for sparse data and index vectors
        static unsigned long get_non_zeros(std::string filename)
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
                unsigned long n_z((unsigned long)atol(n_z_s.c_str()));

                return n_z;

            }
            else
                throw honei::InternalError("Unable to open MatrixMarket file.");
        }

        ///Read in sparse data only
        template<typename DT_>
            static void read_matrix(std::string filename, DenseVector<unsigned long> & row_indices,
                    DenseVector<unsigned long> & column_indices,
                    DenseVector<DT_> & data)
            {
                unsigned long columns, rows, non_data_lines, non_zeros;
                get_sizes(filename, rows, columns, non_zeros, non_data_lines);

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

                    row_indices[i] = r - 1;
                    column_indices[i] = c -1;
                    data[i] = n_z;

                }
                file.close();
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
};
template<>
class MatrixIO<io_formats::ELL>
{
    public:
    static void write_matrix(std::string output, SparseMatrixELL<double> smatrix)
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

    template <typename DT_>
    static SparseMatrixELL<DT_> read_matrix(std::string input, DT_ datatype)
    {
            FILE* file;
            file = fopen(input.c_str(), "rb");
            if (file == NULL)
                throw InternalError("File "+input+" not found!");
            uint64_t size;
            uint64_t rows;
            uint64_t columns;
            uint64_t stride;
            uint64_t num_cols_per_row;
            fread(&size, sizeof(uint64_t), 1, file);
            fread(&rows, sizeof(uint64_t), 1, file);
            fread(&columns, sizeof(uint64_t), 1, file);
            fread(&stride, sizeof(uint64_t), 1, file);
            fread(&num_cols_per_row, sizeof(uint64_t), 1, file);
            uint64_t aj[size];
            DenseVector<double> ax(size);
            fread(aj, sizeof(uint64_t), size, file);
            fread(ax.elements(), sizeof(double), size, file);
            fclose(file);
            DenseVector<DT_> axc(size);
            DenseVector<unsigned long> ajc(size);
            unsigned long crows(rows);
            unsigned long ccolumns(columns);
            unsigned long cstride(stride);
            unsigned long cnum_cols_per_row(num_cols_per_row);
            convert<tags::CPU>(axc, ax);
            for (unsigned long i(0) ; i < size ; ++i)
            {
                ajc[i] = aj[i];
            }
            SparseMatrixELL<DT_> smatrix(crows, ccolumns, cstride, cnum_cols_per_row, ajc, axc);
            return smatrix;
    }
};
#endif
