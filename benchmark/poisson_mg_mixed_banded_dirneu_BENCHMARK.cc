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

//#define SOLVER_VERBOSE

#include <honei/math/multigrid.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <benchmark/benchmark.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/math/endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace std;

template <typename ITag_, typename OTag_, typename DT1_>
class PoissonBenchmarkMGBandedQ1Mixed:
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        PoissonBenchmarkMGBandedQ1Mixed(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
    {
        register_tag(OTag_::name);
        _size  = size;
        _count = count;
    }

        virtual void run()
        {
            unsigned long _root_n(_size);
            unsigned long n(_root_n * _root_n);
            MGInfo<float> info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;
            DenseVector<unsigned long> mask(8);

            info.macro_border_mask = new DenseVector<unsigned long>(8);
            for(unsigned long i(0); i < 8; ++i)
            {
                (*info.macro_border_mask)[i] = 2;
            }
            //set Neumann boundaries:
            (*info.macro_border_mask)[5] =1;

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
            }

            info.n_max_iter = 2;
            info.initial_zero = false;
            info.tolerance = 1e-2;
            info.convergence_check = false;

            info.n_pre_smooth = 2;
            info.n_post_smooth = 2;
            info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.tolerance_coarse = 1e-2;
            info.adapt_correction_factor = 1.;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
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
                unsigned long size = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
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
                unsigned long N = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                DenseVector<float> band(N);
                BandedMatrixQ1<float> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                DenseVector<float> current_rhs(N);


                FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_matrix);

                FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_rhs);

                info.rhs.push_back(current_rhs);
                info.a.push_back(current_matrix);
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
                info.x[i] = null.copy();
            }
            //SET DIAG_INVERTED:
            for (unsigned long i(0) ; i <= info.max_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<float> scaled_diag_inverted(info.a[i].band(DD).copy());
                ElementInverse<OTag_>::value(scaled_diag_inverted);
                Scale<OTag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted[i] = scaled_diag_inverted.copy();
            }

            DenseVector<DT1_> null(info.rhs[info.max_level].size() , DT1_(0));
            BandedMatrixQ1<DT1_> A(info.rhs[info.max_level].size() , null.copy(), null.copy() , null.copy(), null.copy() , null.copy(), null.copy(), null.copy(), null.copy(), null.copy());
            FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(A);
            DenseVector<DT1_> RHS( info.rhs[info.max_level].size(), DT1_(0.));
            FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(RHS);
            DenseVector<DT1_> result(n, DT1_(0));
            for (unsigned long i(0) ; i < 1 ; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 1 ; ++j)
                        {
                        (result = Multigrid<ITag_, OTag_, JAC, CYCLE::V, MIXED >::value(A, RHS, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info));
                        }
                        );
            }
            evaluate();
        }
};
#ifdef HONEI_SSE
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_2_cpu("MG mixed CPU/CPU L2", 5, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_3_cpu("MG mixed CPU/CPU L3", 9, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_4_cpu("MG mixed CPU/CPU L4", 17, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_5_cpu("MG mixed CPU/CPU L5", 33, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_6_cpu("MG mixed CPU/CPU L6", 65, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_7_cpu("MG mixed CPU/CPU L7", 129, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_8_cpu("MG mixed CPU/CPU L8", 257, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_9_cpu("MG mixed CPU/CPU L9", 513, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::CPU::SSE, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_10_cpu("MG mixed CPU/CPU L10", 1025, 1);
#ifdef HONEI_CUDA
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_2("MG mixed GPU/CPU L2", 5, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_3("MG mixed GPU/CPU L3", 9, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_4("MG mixed GPU/CPU L4", 17, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_5("MG mixed GPU/CPU L5", 33, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_6("MG mixed GPU/CPU L6", 65, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_7("MG mixed GPU/CPU L7", 129, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_8("MG mixed GPU/CPU L8", 257, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_9("MG mixed GPU/CPU L9", 513, 1);
PoissonBenchmarkMGBandedQ1Mixed<tags::GPU::CUDA, tags::CPU::SSE, double> poisson_bench_mg_banded_mixed_10("MG mixed GPU/CPU L10", 1025, 1);
#endif
#endif
