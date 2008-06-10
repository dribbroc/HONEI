/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
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

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/lbm/equilibrium_distribution.hh>

using namespace std;
using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_, typename Dir_>
class EquilibriumDistributionBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        EquilibriumDistributionBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> h(1000ul, 1000ul, DataType_(1.23456));
                DenseMatrix<DataType_> u(1000ul, 1000ul, DataType_(1.23456));
                DenseMatrix<DataType_> v(1000ul, 1000ul, DataType_(1.23456));
                DataType_ g(9.81);
                DataType_ e(1.);
                DataType_ e_u(2.);
                DataType_ e_v(2.);

                DenseMatrix<DataType_> result_1(1000ul, 1000ul);
                DenseMatrix<DataType_> result_2(1000ul, 1000ul);
                DenseMatrix<DataType_> result_3(1000ul, 1000ul);

                BENCHMARK((EquilibriumDistribution<Tag_, lbm_applications::LABSWE, Dir_>::value(result_1, h, u, v, g, e)));

            }
            evaluate();
        }
};

EquilibriumDistributionBench<tags::CPU, float, D2Q9::DIR_0> solver_bench_float_0("EquilibriumDistribution DIR_0 Benchmark - size: 1000, float", 1000, 100);
EquilibriumDistributionBench<tags::CPU, double, D2Q9::DIR_0> solver_bench_double_0("EquilibriumDistribution DIR_0 Benchmark - size: 1000, double", 1000, 100);

#ifdef HONEI_SSE
EquilibriumDistributionBench<tags::CPU::SSE, float, D2Q9::DIR_0> solver_bench_float_0_sse("EquilibriumDistribution DIR_0 Benchmark - size: 1000, float SSE", 1000, 100);
EquilibriumDistributionBench<tags::CPU::SSE, double, D2Q9::DIR_0> solver_bench_double_0_sse("EquilibriumDistribution DIR_0 Benchmark - size: 1000, double SSE", 1000, 100);
#endif

template <typename Tag_, typename DataType_, typename Dir_>
class EquilibriumDistributionODDBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        EquilibriumDistributionODDBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> h(1000ul, 1000ul, DataType_(1.23456));
                DenseMatrix<DataType_> u(1000ul, 1000ul, DataType_(1.23456));
                DenseMatrix<DataType_> v(1000ul, 1000ul, DataType_(1.23456));
                DataType_ g(9.81);
                DataType_ e(1.);
                DataType_ e_u(2.);
                DataType_ e_v(2.);

                DenseMatrix<DataType_> result_1(1000ul, 1000ul);
                DenseMatrix<DataType_> result_2(1000ul, 1000ul);
                DenseMatrix<DataType_> result_3(1000ul, 1000ul);

                BENCHMARK((EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_ODD>::value(result_2, h, u, v, g, e, e_u, e_v)));

            }
            evaluate();
        }
};

EquilibriumDistributionODDBench<tags::CPU, float, D2Q9::DIR_ODD> solver_bench_float_1("EquilibriumDistribution DIR_ODD Benchmark - size: 1000, float", 1000, 100);
EquilibriumDistributionODDBench<tags::CPU, double, D2Q9::DIR_ODD> solver_bench_double_1("EquilibriumDistribution DIR_ODD Benchmark - size: 1000, double", 1000, 100);

#ifdef HONEI_SSE
EquilibriumDistributionODDBench<tags::CPU::SSE, float, D2Q9::DIR_ODD> solver_bench_float_1_sse("EquilibriumDistribution DIR_ODD Benchmark - size: 1000, float SSE", 1000, 100);
EquilibriumDistributionODDBench<tags::CPU::SSE, double, D2Q9::DIR_ODD> solver_bench_double_1_sse("EquilibriumDistribution DIR_ODD Benchmark - size: 1000, double SSE", 1000, 100);
#endif


