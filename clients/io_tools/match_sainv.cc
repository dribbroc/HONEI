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
#include <honei/math/sainv.hh>

using namespace honei;

int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        std::cout<<"Usage 'matchsainv match-file ell-file sainv-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if (sizeof(unsigned long) != sizeof(uint64_t))
    {
        std::cout<<"Dont run this on 32bit machines"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in ell file matrix
    std::string match(argv[1]);
    std::string input(argv[2]);
    std::string output(argv[3]);

    SparseMatrixELL<double> smatch(MatrixIO<io_formats::ELL>::read_matrix(match, double(0)));
    unsigned long ue_target = smatch.used_elements();
    SparseMatrixELL<double> smatrix(MatrixIO<io_formats::ELL>::read_matrix(input, double(0)));
    SparseMatrix<double> ssmatrix(smatrix);

    double eps = 0.1;
    unsigned long ue(0);

    SparseMatrix<double> ssainv(1,1,1);
    while (ue < ue_target)
    {
#ifdef HONEI_SSE
        ssainv=(SAINV<tags::CPU::MultiCore::SSE>::value(ssmatrix, eps));
#else
        ssainv=(SAINV<tags::CPU>::value(ssmatrix, eps));
#endif
        ue = ssainv.used_elements();
        eps -= 0.01;
    }

    eps+=0.02;
    ue = 0;

    while (ue < ue_target)
    {
#ifdef HONEI_SSE
        ssainv=(SAINV<tags::CPU::MultiCore::SSE>::value(ssmatrix, eps));
#else
        ssainv=(SAINV<tags::CPU>::value(ssmatrix, eps));
#endif
        ue = ssainv.used_elements();
        eps -= 0.001;
    }

    eps+=0.002;
    ue = 0;

    while (ue < ue_target)
    {
#ifdef HONEI_SSE
        ssainv=(SAINV<tags::CPU::MultiCore::SSE>::value(ssmatrix, eps));
#else
        ssainv=(SAINV<tags::CPU>::value(ssmatrix, eps));
#endif
        ue = ssainv.used_elements();
        eps -= 0.0001;
    }

    std::cout<<"UE Goal: "<<ue_target<<std::endl;
    std::cout<<"Used Elements: "<<ue<<std::endl;
    SparseMatrixELL<double> sainv(ssainv);

    // Write out spai matrix
    MatrixIO<io_formats::ELL>::write_matrix(output, sainv);

    return EXIT_SUCCESS;
}
