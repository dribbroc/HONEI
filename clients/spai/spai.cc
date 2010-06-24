/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 20010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#include <time.h>
#include <honei/math/matrix_io.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/norm.hh>
#include <honei/la/difference.hh>
#include <honei/la/product.hh>

using namespace honei;

int main(int argc, char ** argv)
{
    srand(time(NULL));
    /*if (argc != 3)
    {
        std::cout<<"Usage 'm2ell m-file ell-file'"<<std::endl;
        return EXIT_FAILURE;
    }
    std::string input(argv[1]);*/
    std::string input("../../honei/math/testdata/5pt_10x10.ell");
    std::string spai_input("../../honei/math/testdata/5pt_10x10_spai.ell");

    SparseMatrixELL<double> smatrixell = MatrixIO<io_formats::ELL>::read_matrix(input, double(1));
    SparseMatrix<double> smatrix(smatrixell);
    SparseMatrixELL<double> spai_ref = MatrixIO<io_formats::ELL>::read_matrix(spai_input, double(1));
    SparseMatrix<double> sspai(spai_ref);
    SparseMatrix<double> ident(smatrix.rows(), smatrix.columns(), 1);
    for (unsigned long i(0) ; i < ident.rows() ; ++i)
    {
        ident(i, i) = 1;
    }

    SparseMatrix<double> c(smatrix.rows(), smatrix.columns());

    for (unsigned long i(0) ; i < c.rows() ; ++i)
    {
        c(i, i) = 1. / smatrix(i, i);
    }
    SparseMatrix<double> temp(c.rows(), c.columns());
    double min = Norm<vnt_l_one, false, tags::CPU::SSE>::value(Difference<tags::CPU>::value(temp, ident, Product<tags::CPU>::value(smatrix, c)));
    std::cout<<"C jacobi style "<<min<<std::endl;
    SparseMatrix<double> temp2(c.rows(), c.columns());
    double min2 = Norm<vnt_l_one, false, tags::CPU::SSE>::value(Difference<tags::CPU>::value(temp2, ident, Product<tags::CPU>::value(smatrix, sspai)));
    std::cout<<"SPAI ref "<<min2<<std::endl;


    for (unsigned long iter(0) ; iter < 10 ; ++iter)
    {
        SparseMatrix<double> ttemp(c.rows(), c.columns());
        SparseMatrix<double> tc(c.copy());
        unsigned long changes((rand()%100)+1);
        //unsigned long changes(1);
        for (unsigned long j(0) ; j < changes ; ++j)
        {
            tc(rand()%tc.rows(), rand()%tc.columns()) = double(rand()%tc.size()) / double(rand()%tc.size());
        }
        double tmin = Norm<vnt_l_one, false, tags::CPU::SSE>::value(Difference<tags::CPU>::value(ttemp, ident, Product<tags::CPU>::value(smatrix, tc)));
        if (tmin < min)
        {
            c = tc;
            min = tmin;
            std::cout<<"after "<<iter<<" iters: min = "<<min<<std::endl;
            //break;
        }
    }

    std::cout<<"min nach 100 iters: "<<min<<std::endl;
    //std::cout<<c;

    return EXIT_SUCCESS;
}
