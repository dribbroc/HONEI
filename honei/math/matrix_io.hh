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
#include <string>
#include <honei/la/dense_matrix.hh>

using namespace honei;
class MatrixIO
{
    private:
        static void get_sizes(std::string filename, unsigned long & c,
                                                    unsigned long & r,
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
                            c_s.append(1, line[i]);
                        }
                        line.erase(0, first_blank + 1);

                        std::string::size_type second_blank(line.find_first_of(" "));
                        for(unsigned long i(0) ; i < second_blank ; ++i)
                        {
                            r_s.append(1, line[i]);
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
        template<typename DT_>
            static DenseMatrix<DT_> read_matrix(std::string filename, DT_ base)
            {
                unsigned long columns, rows, non_zeros, non_data_lines;
                get_sizes(filename, columns, rows, non_zeros, non_data_lines);

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
                    for(unsigned long i(0) ; i < first_blank ; ++i)
                    {
                        c_s.append(1, line[i]);
                    }
                    line.erase(0, first_blank + 1);

                    std::string::size_type second_blank(line.find_first_of(" "));
                    for(unsigned long i(0) ; i < second_blank ; ++i)
                    {
                        r_s.append(1, line[i]);
                    }
                    line.erase(0, second_blank + 1);

                    std::string::size_type eol(line.length());
                    for(unsigned long i(0) ; i < eol ; ++i)
                    {
                        n_z_s.append(1, line[i]);
                    }

                    unsigned long c = (unsigned long)atol(c_s.c_str());
                    unsigned long r = (unsigned long)atol(r_s.c_str());
                    DT_ n_z = (DT_)atof(n_z_s.c_str());

                    result[r - 1][c - 1] = n_z;
                }

                return result;
            }
};
#endif
