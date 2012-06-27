/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
    if (argc != 5)
    {
        std::cout<<"Usage 'ellexpand file target-file new-rows new-columns'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in m file matrix
    std::string input(argv[1]);
    std::string output(argv[2]);
    unsigned long row(atoi(argv[3]));
    unsigned long column(atoi(argv[4]));

    SparseMatrixELL<double> tsmatrix(MatrixIO<io_formats::ELL>::read_matrix(input, double(0)));

    SparseMatrix<double> temp(tsmatrix);
    SparseMatrix<double> ssmatrix(row, column);
    for (SparseMatrix<double>::NonZeroConstElementIterator iter(temp.begin_non_zero_elements()) ; iter != temp.end_non_zero_elements() ; ++iter)
    {
        ssmatrix(iter.row(), iter.column(), *iter);
    }
    SparseMatrixELL<double> smatrix(ssmatrix);
    MatrixIO<io_formats::ELL>::write_matrix(output, smatrix);

    return EXIT_SUCCESS;
}
