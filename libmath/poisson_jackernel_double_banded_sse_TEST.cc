/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
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

#include <libmath/jacobi_kernel.hh>
#include <libmath/conjugate_gradients.hh>
#include <libmath/iterative_refinement.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class PoissonTestJACKernelBandedDoubleSSE:
    public BaseTest
{
    public:
        PoissonTestJACKernelBandedDoubleSSE(const std::string & tag) :
            BaseTest("Poisson test for itrerative LES solvers , JAC Kernel(banded system) <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            int n;
            FILE* file;

            double* dd;

            double* ll;
            double* ld;
            double* lu;
            double* dl;
            double* du;
            double* ul;
            double* ud;
            double* uu;

            double* b;
            double* ana_sol;
            double* ref_sol;

            file = fopen("ehq.1.1.1.1.bin", "rb");
            fread(&n, sizeof(int), 1, file);
            dd = new double[n];
            ll = new double[n];
            ld = new double[n];
            lu = new double[n];
            dl = new double[n];
            du = new double[n];
            ul = new double[n];
            ud = new double[n];
            uu = new double[n];
            b = new double[n];
            ana_sol = new double[n];
            ref_sol = new double[n];

            fread(dd, sizeof(double), n, file);
            fread(ll, sizeof(double), n, file);
            fread(ld, sizeof(double), n, file);
            fread(lu, sizeof(double), n, file);
            fread(dl, sizeof(double), n, file);
            fread(du, sizeof(double), n, file);
            fread(ul, sizeof(double), n, file);
            fread(ud, sizeof(double), n, file);
            fread(uu, sizeof(double), n, file);
            fread(b,  sizeof(double), n, file);
            fread(ana_sol, sizeof(double), n, file);
            fread(ref_sol, sizeof(double), n, file);
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
            DenseVector<double> b_v(n, double(0));
            DenseVector<double> ana_sol_v(n, double(0));
            DenseVector<double> ref_sol_v(n, double(0));
            for(unsigned long i = 0; i < n; ++i)
            {
                dd_v[i] = dd[i];
                ll_v[i] = ll[i];
                ld_v[i] = ld[i];
                lu_v[i] = lu[i];
                dl_v[i] = dl[i];
                du_v[i] = du[i];
                ul_v[i] = ul[i];
                ud_v[i] = ud[i];
                uu_v[i] = uu[i];
                b_v[i] = b[i];
                ana_sol_v[i] = ana_sol[i];
                ref_sol_v[i] = ref_sol[i];
            }
            //std::cout<<dd[4]<<endl;
            //std::cout<<dd_v<<endl;


            unsigned int root_n = (unsigned int)sqrt(n);
            BandedMatrix<double> A(n,dd_v.copy());
            //std::cout<<A.band(0)<<endl;
            //A->insert_band(0, dd_v.copy());
            A.insert_band(1, du_v);
            A.insert_band(-1, dl_v);
            A.insert_band(root_n, ud_v);
            A.insert_band(root_n+1, uu_v);
            A.insert_band(root_n-1, ul_v);
            A.insert_band(-root_n, ld_v);
            A.insert_band(-root_n-1, ll_v );
            A.insert_band(-root_n+1, lu_v);
            //New:
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(ref_sol_v);
            DenseVector<DT1_> x(b_v.size(), DT1_(0));
            DenseVector<DT1_> x_last(x.copy());

            DT1_ norm_x_last = DT1_(0);
            DT1_ norm_x = DT1_(1);
            DenseVector<DT1_> diag(b_v.size(), DT1_(0));

            DenseVector<DT1_> diag_inverted(b_v.size(), DT1_(0));

            BandedMatrix<DT1_> difference(A.copy());
            ///Create Diagonal, invert, compute difference on the fly.
            for(unsigned long i =0; i < diag.size(); ++i)
            {
                diag[i] = A.band(0)[i];
                if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                {
                    diag_inverted[i] = DT1_(1) / diag[i];
                }
                else
                {
                    diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                }
            }

            DenseVector<DT1_> zeros(b_v.size(), DT1_(0));
            difference.insert_band(0, zeros);
            Scale<tags::CPU>::value(difference, DT1_(-1));

            DT1_ konv_rad = std::numeric_limits<DT1_>::epsilon();
            while(fabs(norm_x - norm_x_last) > konv_rad)
            {

                JacobiKernel<Tag_>::value(b_v, x, diag_inverted, difference);
                norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                x_last = x.copy();
            }
            cout << "x: " << x << endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, norm_x , 10e-06);

        }
};
#ifdef HONEI_SSE
PoissonTestJACKernelBandedDoubleSSE<tags::CPU::SSE, double> poisson_test_jac_banded_double_sse("SSE double");
#endif