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

#ifndef MATH_GUARD_VECTOR_IO_HH
#define MATH_GUARD_VECTOR_IO_HH 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <honei/la/dense_vector.hh>
#include <honei/la/algorithm.hh>

using namespace honei;

namespace io_formats
{
    class EXP;
    class DV;
}


template<typename IOFormat_>
class VectorIO
{
};

template<>
class VectorIO<io_formats::EXP>
{
    public:
        static void get_size(std::string filename, unsigned long & n_z, unsigned long & data_begin)
        {
            n_z = 0;
            data_begin = 0;

            std::ifstream file(filename.c_str());
            if(file.is_open())
            {
                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);

                    if(line.find("#", 0) < line.npos)
                    {
                        ++data_begin;
                    }
                    else
                    {
                        ++n_z;
                    }

                }

                file.close();
                --n_z;
            }
            else
                throw honei::InternalError("Unable to open Vector file.");
        }

        template<typename DT_>
            static void read_vector(std::string filename, DenseVector<DT_> & data)
            {
                unsigned long non_data_lines, non_zeros;
                get_size(filename, non_zeros, non_data_lines);

                std::ifstream file(filename.c_str());
                file.is_open();

                for(unsigned long i(0) ; i < non_data_lines ; ++i)
                {
                    std::string line;
                    std::getline(file, line);
                }

                for(unsigned long j(0) ; j < non_zeros ; ++j)
                {
                    std::string line;
                    std::getline(file, line);

                    std::string n_z_s;

                    std::string::size_type first_digit(line.find_first_not_of(" "));
                    line.erase(0, first_digit);
                    std::string::size_type eol(line.length());
                    for(unsigned long i(0) ; i < eol ; ++i)
                    {
                        n_z_s.append(1, line[i]);
                    }

                    DT_ n_z = (DT_)atof(n_z_s.c_str());

                    data[j] = n_z;

                }
                file.close();

            }
};

template<>
class VectorIO<io_formats::DV>
{
    public:
    static void write_vector(std::string output, DenseVector<double> dv)
    {
            FILE* file;
            file = fopen(output.c_str(), "wb");
            uint64_t size(dv.size());
            fwrite(&size, sizeof(uint64_t), 1, file);
            fwrite(dv.elements(), sizeof(double), size, file);
            fclose(file);
    }

    template <typename DT_>
    static DenseVector<DT_> read_vector(std::string input, DT_ datatype)
    {
            FILE* file(NULL);
            file = fopen(input.c_str(), "rb");
            if (file == NULL)
                throw InternalError("File "+input+" not found!");
            uint64_t size;
            fread(&size, sizeof(uint64_t), 1, file);
            DenseVector<double> ax(size);
            fread(ax.elements(), sizeof(double), size, file);
            fclose(file);
            DenseVector<DT_> axc(size);
            convert<tags::CPU>(axc, ax);
            return axc;
    }

};
#endif
