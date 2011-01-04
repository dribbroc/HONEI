/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

#include <string>
#endif

#include <honei/math/iterative_refinement.hh>
#include <honei/math/endian_swap.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>

class PoissonIteRefPCGBenchDouble :
    public Benchmark
{
    public:
        PoissonIteRefPCGBenchDouble(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
        }

        virtual void run()
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

            long root_n = (long)sqrt(n);
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

            //std::cout<<A.band(0)[0] * double(1) << endl;

            //std::cout<< n << " " << A << " "<< root_n<<endl;
            DenseVector<double> result(n, double(0));
            BENCHMARK((IterativeRefinement<methods::PCG::JAC, Tag_>::value(A, b_v, result, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon())));
            evaluate();
        }
};
PoissonIteRefPCGBenchDouble<tags::CPU, double> poisson_irpcg_bench_double("Poisson IteRefPCG benchmark double CPU");
#ifdef HONEI_SSE
PoissonIteRefPCGBenchDouble<tags::CPU::SSE, double> poisson_irpcg_bench_double_sse("Poisson IteRefPCG benchmark double SSE");
#endif
#ifdef HONEI_CELL
PoissonIteRefPCGBenchDouble<tags::Cell, double> poisson_irpcg_bench_double_cell("Poisson IteRefPCG benchmark double Cell");
#endif
