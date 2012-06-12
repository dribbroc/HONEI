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
#ifndef MPI_GUARD_VECTOR_IO_MPI_HH
#define MPI_GUARD_VECTOR_IO_MPI_HH 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <honei/math/vector_io.hh>
#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/algorithm.hh>

using namespace honei;

template<typename IOFormat_>
class VectorIOMPI
{
};

template<>
class VectorIOMPI<io_formats::EXP>
{
    public:
        template<typename DT_>
            static DenseVectorMPI<DT_> read_vector(std::string filename, DT_)
            {
                unsigned long _rank(mpi::mpi_comm_rank());
                if (_rank == 0)
                {
                    unsigned long global_size(0);

                    std::ifstream file(filename.c_str());
                    if (! file.is_open())
                        throw honei::InternalError("Unable to open Vector file " + filename);
                    while(!file.eof())
                    {
                        std::string line;
                        std::getline(file, line);
                        if(line.find("#", 0) < line.npos)
                            continue;
                        if(file.eof())
                            break;

                        ++global_size;
                    }
                    mpi::mpi_bcast(&global_size, 1, 0);


                    DenseVector<DT_> data(global_size);
                    unsigned long index(0);

                    file.clear();
                    file.seekg (0, std::ios::beg);

                    while(!file.eof())
                    {
                        std::string line;
                        std::getline(file, line);
                        if(line.find("#", 0) < line.npos)
                            continue;
                        if(file.eof())
                            break;


                        std::string n_z_s;

                        std::string::size_type first_digit(line.find_first_not_of(" "));
                        line.erase(0, first_digit);
                        std::string::size_type eol(line.length());
                        for(unsigned long i(0) ; i < eol ; ++i)
                        {
                            n_z_s.append(1, line[i]);
                        }

                        DT_ n_z = (DT_)atof(n_z_s.c_str());

                        data[index] = n_z;

                        ++index;
                    }
                    file.close();

                    mpi::mpi_bcast(data.elements(), data.size(), 0);
                    DenseVectorMPI<DT_> result(data);
                    return result;
                }
                else
                {
                    unsigned long global_size;
                    mpi::mpi_bcast(&global_size, 1, 0);

                    DenseVector<DT_> data(global_size);
                    mpi::mpi_bcast(data.elements(), data.size(), 0);
                    DenseVectorMPI<DT_> result(data);
                    return result;
                }

            }
};
#endif
