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
#include <honei/math/matrix_io.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/banded_matrix_qx.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <cmath>
#include <fstream>

using namespace honei;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'feast2ell feast-bin-file ell-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in feast-bin file matrix
    int n;
    std::string input(argv[1]);
    std::string output(argv[2]);

    FILE* file = fopen(input.c_str(), "rb");
    if (1 != fread(&n, sizeof(int), 1, file))
        throw InternalError("IO Error!");

    double * dd = new double[n];
    double * ll = new double[n];
    double * ld = new double[n];
    double * lu = new double[n];
    double * dl = new double[n];
    double * du = new double[n];
    double * ul = new double[n];
    double * ud = new double[n];
    double * uu = new double[n];
    //double * b = new double[n];
    //double * ana_sol = new double[n];
    //double * ref_sol = new double[n];

    if (n != (int)fread(dd, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ll, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ld, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(lu, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(dl, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(du, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ul, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ud, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(uu, sizeof(double), n, file))
        throw InternalError("IO Error!");
    /*if (n != (int)fread(b,  sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ana_sol, sizeof(double), n, file))
        throw InternalError("IO Error!");
    if (n != (int)fread(ref_sol, sizeof(double), n, file))
        throw InternalError("IO Error!");*/
    fclose(file);

    DenseVector<double> dd_v(n, double(0));
    DenseVector<double> ll_v(n, double(0));
    DenseVector<double> ld_v(n, double(0));
    DenseVector<double> lu_v(n, double(0));
    DenseVector<double> dl_v(n, double(0));
    DenseVector<double> du_v(n, double(0));
    DenseVector<double> ul_v(n, double(0));
    DenseVector<double> ud_v(n, double(0));
    DenseVector<double> uu_v(n, double(0));
    //DenseVector<double> b_v(n, double(0));
    //DenseVector<double> ana_sol_v(n, double(0));
    //DenseVector<double> ref_sol_v(n, double(0));
    for(int i = 0; i < n; ++i)
    {
        dd_v[i] = double(dd[i]);
        ll_v[i] = ll[i];
        ld_v[i] = ld[i];
        lu_v[i] = lu[i];
        dl_v[i] = dl[i];
        du_v[i] = du[i];
        ul_v[i] = ul[i];
        ud_v[i] = ud[i];
        uu_v[i] = uu[i];
        //b_v[i] = b[i];
        //ana_sol_v[i] = ana_sol[i];
        //ref_sol_v[i] = ref_sol[i];
    }

    long root_n = (long)sqrt(n);
    BandedMatrix<double> A(n,dd_v.copy());
    A.insert_band(1, du_v);
    A.insert_band(-1, dl_v);
    A.insert_band(root_n, ud_v);
    A.insert_band(root_n+1, uu_v);
    A.insert_band(root_n-1, ul_v);
    A.insert_band(-root_n, ld_v);
    A.insert_band(-root_n-1, ll_v );
    A.insert_band(-root_n+1, lu_v);

    BandedMatrixQx<Q1Type, double> Aq1(A);
    SparseMatrix<double> sm(Aq1);
    SparseMatrixELL<double> smell(sm);

    // Write out ell file matrix
    MatrixIO<io_formats::ELL>::write_matrix(output, smell);

    return EXIT_SUCCESS;
}
