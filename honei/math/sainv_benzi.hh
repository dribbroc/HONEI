/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#pragma once
#ifndef LIBMATH_GUARD_SAINV_BENZI_HH
#define LIBMATH_GUARD_SAINV_BENZI_HH 1

#include <honei/util/tags.hh>
#include <honei/la/sparse_matrix.hh>


// Based on "Robust Approximate Inverse Preconditioning for the Conjugate Gradients Method" by Benzi et al.
// http://www2.cs.cas.cz/~tuma/sparslab.html

extern "C" void ainvsr2_(int & msglvl, int & msgunit, int & n, int * ia, int * ja, double * a, int * ip, int * jp, double * ap,
     int & size_p, int & size_c, int & size_r, double & diagtol,
     double & drfl, double & mi, int & diag_one, int & droptyp, int & imodif, int & fill, int & fillmax,
     int & ifillmax, int & garrow, int & garcol, int & info);

namespace honei
{
    struct SAINV_BENZI
    {
        static SparseMatrix<double> value(const SparseMatrix<double> & A, double tolerance = 1e-1)
        {
            int size_p=1000000;
            int size_r=1000000;
            int size_c=1000000;
            double diagtol=1.1e-16;
            int droptyp=0;
            double mi=0.1;
            int diag_one=1; //1==scale diag by pivot
            int msglvl=0; //3 for verbosive output
            int msgunit=0;

            int n=A.rows();
            int  * ia = new int[n + 1];
            int * ja = new int[A.used_elements()];
            double * aa = new double[A.used_elements()];
            int ue(1);
            int ij(0);
            for (unsigned long i(0) ; i < A.rows() ; ++i)
            {
                ia[i] = ue;
                for (unsigned long j(0) ; j < A[i].used_elements() ; ++j)
                {
                    ja[ij] = A[i].indices()[j] + 1;
                    aa[ij] = A[i].elements()[j];
                    ++ij;
                }
                ue += A[i].used_elements();
            }
            ia[n] = ue;


            /*int n=4;
            int  * ia = new int[n + 1];
            int * ja = new int[12];
            double * aa = new double[12];
            ia[0]=1;
            ia[1]=4;
            ia[2]=7;
            ia[3]=10;
            ia[4]=13;
            ja[0]=1;
            ja[1]=2;
            ja[2]=3;
            ja[3]=1;
            ja[4]=2;
            ja[5]=4;
            ja[6]=1;
            ja[7]=3;
            ja[8]=4;
            ja[9]=2;
            ja[10]=3;
            ja[11]=4;
            aa[0]=4;
            aa[1]=-1;
            aa[2]=-1;
            aa[3]=-1;
            aa[4]=4;
            aa[5]=-1;
            aa[6]=-1;
            aa[7]=4;
            aa[8]=-1;
            aa[9]=-1;
            aa[10]=-1;
            aa[11]=4;
            int ue = 13;*/



            int  * ip = new int[ue*2];
            int * jp = new int[ue*100];
            double * ap = new double[ue*100];

            int imodif = 0;
            int fill = 0;
            int fillmax = 0;
            int ifillmax = 0;
            int garrow = 0;
            int garcol = 0;
            int info = 0;


            ainvsr2_(msglvl,msgunit,n,
                    ia,ja,aa,ip,jp,ap,
                    size_p,size_c,size_r,diagtol,
                    tolerance,mi,diag_one,droptyp,imodif,fill,fillmax,
                    ifillmax,garrow,garcol,info);

            //SparseMatrix<double> result(A.rows(), A.columns());
            SparseMatrix<double> result(4, 4);
            for (int i(0) ; i < n ; ++i)
            {
                for (int j(ip[i]-1) ; j < ip[i+1]-1 ; ++j)
                {
                    result(i, jp[j] - 1, ap[j]);
                }
            }
            std::cout<<result;

            delete[] ia;
            delete[] aa;
            delete[] ja;
            delete[] ip;
            delete[] jp;
            delete[] ap;

            return result;
        }
    };
}
#endif
