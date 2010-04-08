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

//#define SOLVER_VERBOSE 1
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

template <typename Tag_, typename DT1_>
class PoissonBenchmarkMGBandedQ1Fixed:
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        PoissonBenchmarkMGBandedQ1Fixed(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
    {
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
    }

        virtual void run()
        {
            unsigned long _root_n(_size);
            unsigned long n(_root_n * _root_n);
            MGInfo<DT1_, BandedMatrixQ1<DT1_> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;
            DenseVector<unsigned long> mask(8);

            info.macro_border_mask = new DenseVector<unsigned long>(8);
            for(int i(0); i < 8; ++i)
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
                default:
                    throw InternalError("Uknown size!");
                    break;
            }

            info.n_max_iter = 16;
            info.initial_zero = false;
            info.tolerance = 1e-8;
            info.convergence_check = true;

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

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                BandedMatrixQ1<DT1_> ac_a(size, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band);
                info.a.push_back(ac_a);
                // iteration vectors
                DenseVector<DT1_> ac_c(size, DT1_(0));
                info.c.push_back(ac_c);
                DenseVector<DT1_> ac_d(size, DT1_(0));
                info.d.push_back(ac_d);
                DenseVector<DT1_> ac_rhs(size, DT1_(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<DT1_> ac_x(size, DT1_(0));
                info.x.push_back(ac_x);

                info.diags_inverted.push_back(dummy_band.copy());
            }


            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                // iteration vectors
                DenseVector<DT1_> ac_c(size, DT1_(0));
                info.c.push_back(ac_c);
                DenseVector<DT1_> ac_d(size, DT1_(0));
                info.d.push_back(ac_d);
                DenseVector<DT1_> ac_x(size, DT1_(0));
                info.x.push_back(ac_x);

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                info.diags_inverted.push_back(dummy_band.copy());
            }

            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                DenseVector<DT1_> band(N);
                BandedMatrixQ1<DT1_> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                DenseVector<DT1_> current_rhs(N);


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

                DenseVector<DT1_> null(size , DT1_(0));
                info.x[i] = null.copy();
            }
            //SET DIAG_INVERTED:
            for (unsigned long i(0) ; i <= info.max_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<DT1_> scaled_diag_inverted(info.a[i].band(DD).copy());
                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted[i] = scaled_diag_inverted.copy();
            }

            DenseVector<DT1_> result(n, DT1_(0));
            for (int i(0) ; i < 1 ; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 1 ; ++j)
                        {
                        (Multigrid<Tag_, Tag_, JAC, CYCLE::V, FIXED >::value(info.a[info.max_level], info.rhs[info.max_level], result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info));
                        }
                        );
            }
            evaluate();
        }
};
#ifdef HONEI_SSE
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_2("MG float SSE L2", 5, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_3("MG float SSE L3", 9, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_4("MG float SSE L4", 17, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_5("MG float SSE L5", 33, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_6("MG float SSE L6", 65, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_7("MG float SSE L7", 129, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_8("MG float SSE L8", 257, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_9("MG float SSE L9", 513, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, float> poisson_bench_mg_banded_float_10("MG float SSE L10", 1025, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_2("MG double SSE L2", 5, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_3("MG double SSE L3", 9, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_4("MG double SSE L4", 17, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_5("MG double SSE L5", 33, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_6("MG double SSE L6", 65, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_7("MG double SSE L7", 129, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_8("MG double SSE L8", 257, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_9("MG double SSE L9", 513, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::CPU::SSE, double> poisson_bench_mg_banded_double_10("MG double SSE L10", 1025, 1);
#endif

#ifdef HONEI_CUDA
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a2("MG float CUDA L2", 5, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a3("MG float CUDA L3", 9, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a4("MG float CUDA L4", 17, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a5("MG float CUDA L5", 33, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a6("MG float CUDA L6", 65, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a7("MG float CUDA L7", 129, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a8("MG float CUDA L8", 257, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a9("MG float CUDA L9", 513, 1);
PoissonBenchmarkMGBandedQ1Fixed<tags::GPU::CUDA, float> poisson_bench_mg_banded_float_a10("MG float CUDA L10", 1025, 1);
#endif
