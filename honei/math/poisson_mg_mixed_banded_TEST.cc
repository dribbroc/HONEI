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

#include <honei/math/multigrid.hh>
#include <honei/math/conjugate_gradients.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename OuterTag_>
class PoissonTestMGBandedQ1Mixed:
    public BaseTest
{
    public:
        PoissonTestMGBandedQ1Mixed(const std::string & tag) :
            BaseTest("Poisson test for itrerative LES solvers , MG (Q1 system) MIXED <" + tag + ">")
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

            std::string file_name(HONEI_SOURCEDIR);
            file_name += "/honei/math/testdata/4225.bin";
            file = fopen(file_name.c_str(), "rb");
            if (1 != (int)fread(&n, sizeof(int), 1, file))
                throw InternalError("IO Error!");

#ifdef HONEI_CELL
            unsigned char b1, b2, b3, b4;
            b1 = n & 255;
            b2 = ( n >> 8 ) & 255;
            b3 = ( n>>16 ) & 255;
            b4 = ( n>>24 ) & 255;
            n = ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
#endif
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
            if (n != (int)fread(b,  sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ana_sol, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ref_sol, sizeof(double), n, file))
                throw InternalError("IO Error!");
            fclose(file);

#ifdef HONEI_CELL
            for(int i(0); i < n; ++i)
            {
                dd[i] = DoubleSwap(dd[i]);
                ll[i] = DoubleSwap(ll[i]);
                ld[i] = DoubleSwap(ld[i]);
                lu[i] = DoubleSwap(lu[i]);
                dl[i] = DoubleSwap(dl[i]);
                du[i] = DoubleSwap(du[i]);
                ul[i] = DoubleSwap(ul[i]);
                ud[i] = DoubleSwap(ud[i]);
                uu[i] = DoubleSwap(uu[i]);
                b[i] = DoubleSwap(b[i]);
                ana_sol[i] = DoubleSwap(ana_sol[i]);
                ref_sol[i] = DoubleSwap(ref_sol[i]);

            }
#endif
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
            for(int i = 0; i < n; ++i)
            {
                dd_v[i] = (double)dd[i];
                ll_v[i] = (double)ll[i];
                ld_v[i] = (double)ld[i];
                lu_v[i] = (double)lu[i];
                dl_v[i] = (double)dl[i];
                du_v[i] = (double)du[i];
                ul_v[i] = (double)ul[i];
                ud_v[i] = (double)ud[i];
                uu_v[i] = (double)uu[i];
                b_v[i] = (double)b[i];
                ana_sol_v[i] = (double)ana_sol[i];
                ref_sol_v[i] = (double)ref_sol[i];
            }
            //std::cout<<dd[4]<<endl;
            //std::cout<<dd_v<<endl;


            BandedMatrixQ1<double> A(n, ll_v.copy(), ld_v.copy(), lu_v.copy(), dl_v.copy(), dd_v.copy(), du_v.copy(), ul_v.copy(), ud_v.copy(), uu_v.copy());
            //std::cout<<A.band(0)<<endl;
            //A->insert_band(0, dd_v.copy());
            //std::cout<<A.band(0)[0] * double(1) << endl;
            //std::cout<< n << " " << A << " "<< root_n<<endl;

            //----Load in POISSON - specific data:----------------------
            MGInfo<float, BandedMatrixQ1<float> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;
            DenseVector<unsigned long> mask(8);

            info.macro_border_mask = new DenseVector<unsigned long>(8);
            for(int i(0); i < 8; ++i)
            {
                (*info.macro_border_mask)[i] = 2;
            }

            info.min_level = 1;
            switch(n)
            {
                case 1050625:
                    {
                        info.max_level = 10;
                    }
                    break;
                case 263169:
                    {
                        info.max_level = 9;
                    }
                    break;
                case 66049:
                    {
                        info.max_level = 8;
                    }
                    break;
                case 16641:
                    {
                        info.max_level = 7;
                    }
                    break;
                case 4225:
                    {
                        info.max_level = 6;
                    }
                    break;
                case 1089:
                    {
                        info.max_level = 5;
                    }
                    break;
                case 289:
                    {
                        info.max_level = 4;
                    }
                    break;
                case 81:
                    {
                        info.max_level = 3;
                    }
                    break;
                case 25:
                    {
                        info.max_level = 2;
                    }
                    break;
                case 9:
                    {
                        info.max_level = 1;
                    }
                    break;
                default:
                    throw InternalError("Uknown size!");
                    break;
            }

            info.n_max_iter = 2;
            info.initial_zero = true;
            info.tolerance = 1e-8;
            info.convergence_check = false;

            info.n_pre_smooth = 2;
            info.n_post_smooth = 2;
            info.n_max_iter_coarse = ((unsigned long)sqrt((float)(pow((float)2 , (float)info.max_level) + 1)*(pow((float)2 , (float)info.max_level) + 1)));
            info.tolerance_coarse = 1e-2;
            info.adapt_correction_factor = 1.;

            //push back dummy matrices/vectors in order not to disturb std::vectors index range:
            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<float> dummy_band(size, float(0));
                BandedMatrixQ1<float> ac_a(size, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band);
                info.a.push_back(ac_a);
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_rhs(size, float(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);
                info.diags_inverted.push_back(dummy_band.copy());
            }
            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size = (unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1));
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);
                DenseVector<float> dummy_band(size, float(0));
                info.diags_inverted.push_back(dummy_band.copy());
            }

            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N = (unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1));
                DenseVector<float> LL_v_2(N);
                DenseVector<float> LD_v_2(N);
                DenseVector<float> LU_v_2(N);
                DenseVector<float> DL_v_2(N);
                DenseVector<float> DD_v_2(N);
                DenseVector<float> DU_v_2(N);
                DenseVector<float> UL_v_2(N);
                DenseVector<float> UD_v_2(N);
                DenseVector<float> UU_v_2(N);
                BandedMatrixQ1<float> current_matrix(N,LL_v_2, LD_v_2, LU_v_2, DL_v_2, DD_v_2, DU_v_2, UL_v_2, UD_v_2, UU_v_2);

                DenseVector<float> current_rhs(N);
                int n_2;

                FILE* file_2;

                double* dd_2;

                double* ll_2;
                double* ld_2;
                double* lu_2;
                double* dl_2;
                double* du_2;
                double* ul_2;
                double* ud_2;
                double* uu_2;
                double* b_2;
                std::string file_path(HONEI_SOURCEDIR);
                file_path += "/honei/math/testdata/" + stringify(DD_v_2.size()) +".bin";
                file_2 = fopen(file_path.c_str(), "rb");
                if (1 != (int)fread(&n_2, sizeof(int), 1, file_2))
                    throw InternalError("IO Error!");

#ifdef HONEI_CELL
                unsigned char b1, b2, b3, b4;
                b1 = n_2 & 255;
                b2 = ( n_2 >> 8 ) & 255;
                b3 = ( n_2>>16 ) & 255;
                b4 = ( n_2>>24 ) & 255;
                n_2 = ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
#endif
                dd_2 = new double[n_2];
                ll_2 = new double[n_2];
                ld_2 = new double[n_2];
                lu_2 = new double[n_2];
                dl_2 = new double[n_2];
                du_2 = new double[n_2];
                ul_2 = new double[n_2];
                ud_2 = new double[n_2];
                uu_2 = new double[n_2];

                b_2 = new double[n];
                if (n_2 != (int)fread(dd_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(ll_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(ld_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(lu_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(dl_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(du_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(ul_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(ud_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(uu_2, sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                if (n_2 != (int)fread(b_2,  sizeof(double), n_2, file_2))
                    throw InternalError("IO Error!");
                fclose(file_2);

#ifdef HONEI_CELL
                for(unsigned long j(0); j < n_2; ++j)
                {
                    dd_2[j] = DoubleSwap(dd_2[j]);
                    ll_2[j] = DoubleSwap(ll_2[j]);
                    ld_2[j] = DoubleSwap(ld_2[j]);
                    lu_2[j] = DoubleSwap(lu_2[j]);
                    dl_2[j] = DoubleSwap(dl_2[j]);
                    du_2[j] = DoubleSwap(du_2[j]);
                    ul_2[j] = DoubleSwap(ul_2[j]);
                    ud_2[j] = DoubleSwap(ud_2[j]);
                    uu_2[j] = DoubleSwap(uu_2[j]);

                    b_2[j] = DoubleSwap(b_2[j]);
                }
#endif
                for(unsigned long j(0); j < DD_v_2.size(); ++j)
                {
                    LL_v_2[j] = (float)ll_2[j];
                    LD_v_2[j] = (float)ld_2[j];
                    LU_v_2[j] = (float)lu_2[j];
                    DL_v_2[j] = (float)dl_2[j];
                    DD_v_2[j] = (float)dd_2[j];
                    DU_v_2[j] = (float)du_2[j];
                    UL_v_2[j] = (float)ul_2[j];
                    UD_v_2[j] = (float)ud_2[j];
                    UU_v_2[j] = (float)uu_2[j];
                    current_rhs[j] = (float)b_2[j];
                }
                info.rhs.push_back(current_rhs);
                info.a.push_back(current_matrix);

            }

            //clear rhs data on lower than max_level
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1)));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
                info.x[i] = null.copy();
            }
//SET DIAG_INVERTED:
            for (unsigned long i(0) ; i <= info.max_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((double)2, (double)i) + 1) * ((unsigned long)pow((double)2, (double)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<float> scaled_diag_inverted(info.a[i].band(DD).copy());
                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted[i] = scaled_diag_inverted.copy();
            }
            //--------End loading of data----------------------------------
            DenseVector<double> result(n, double(0));
            Multigrid<Tag_, OuterTag_, methods::NONE, methods::JAC, methods::CYCLE::V, methods::MIXED >::value(A, b_v, result, (unsigned long)11, std::numeric_limits<double>::epsilon(), info);

            result.lock(lm_read_only);
            for(int i = 0; i < n; i++)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(ref_sol[i], result[i], 1e-04);
            }
            result.unlock(lm_read_only);
            DenseVector<double> x(n, double(0));
            Difference<OuterTag_>::value(result, ana_sol_v);
            Difference<OuterTag_>::value(x, result);
            double norm = Norm<vnt_l_two, false, OuterTag_>::value(x);
            std::cout << "L2: " << norm << std::endl;
            //TEST_CHECK(true);
        }
};
PoissonTestMGBandedQ1Mixed<tags::CPU, tags::CPU> poisson_test_mg_bandedMixed("Mixed , CPU/CPU");
#ifdef HONEI_SSE
PoissonTestMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE> sse_poisson_test_mg_bandedMixed("Mixed, SSE/SSE");
#endif
#ifdef HONEI_CELL
//PoissonTestMGBandedQ1Mixed<tags::Cell, tags::Cell> cell_poisson_test_mg_bandedMixed("Mixed CELL/CELL");
#endif
#if defined HONEI_CUDA && defined HONEI_SSE
PoissonTestMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE> cuda_sse_poisson_test_mg_bandedMixed("Mixed CUDA/SSE");
#endif
#ifdef HONEI_CUDA
PoissonTestMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU> cuda_cpu_poisson_test_mg_bandedMixed("Mixed CUDA/CPU");
#endif
#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
PoissonTestMGBandedQ1Mixed<tags::GPU::CUDA, tags::GPU::CUDA> cuda_cuda_poisson_test_mg_bandedMixed("Mixed CUDA/CUDA");
#endif
