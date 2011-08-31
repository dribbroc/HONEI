/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of HONEI. HONEI is free software;
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


#include <iostream>
#include <honei/math/vector_io.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <fstream>

using namespace honei;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'feast2dv feast-bin-file dv-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in feast-bin file matrix
    int n;
    std::string input(argv[1]);
    std::string output(argv[2]);

    FILE* file = fopen(input.c_str(), "rb");
    if (1 != fread(&n, sizeof(int), 1, file))
        throw InternalError("IO Error!");
    double * b = new double[n];

    if (n != (int)fread(b,  sizeof(double), n, file))
        throw InternalError("IO Error!");
    fclose(file);

    DenseVector<double> b_v(n, double(0));
    for(int i = 0; i < n; ++i)
    {
        b_v[i] = b[i];
    }


    // Write out exp vector file
    VectorIO<io_formats::DV>::write_vector(output, b_v);

    return EXIT_SUCCESS;
}
