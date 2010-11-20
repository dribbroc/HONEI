/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#include <honei/math/matrix_io.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>

using namespace honei;

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cout<<"Usage 'ellinfo mtx-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in m file matrix
    std::string input(argv[1]);
    SparseMatrix<double> tsmatrix(MatrixIO<io_formats::ELL>::read_matrix(input, double(0)));

    std::cout<<"ELL Matrix info for " + input << std::endl;
    std::cout<<"Rows: " << tsmatrix.rows() << std::endl;
    std::cout<<"Columns: " << tsmatrix.columns() << std::endl;
    std::cout<<"Non Zero Elements: " << tsmatrix.used_elements() << std::endl;
    return EXIT_SUCCESS;
}
