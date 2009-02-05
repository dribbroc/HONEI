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

#include <honei/lbm/source.hh>

using namespace std;
using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class SourceBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SourceBench(const std::string & id, unsigned long size, int count) :
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
                DenseMatrix<DataType_> db(1000ul, 1000ul, DataType_(0.));
                DataType_ g(9.81);
                DenseMatrix<DataType_> result(1000ul, 1000ul);

                BENCHMARK((Source<Tag_, lbm_applications::LABSWE, lbm_force::SIMPLE, lbm_source_schemes::BASIC>::value(result, h, db, g)));

            }
            evaluate();
        }
};

SourceBench<tags::CPU, float> solver_bench_float_0("Source   Benchmark - size: 1000, float", 1000, 100);
SourceBench<tags::CPU, double> solver_bench_double_0("Source   Benchmark - size: 1000, double", 1000, 100);
#ifdef HONEI_SSE
SourceBench<tags::CPU::SSE, float > solver_bench_float_0_sse("Source   Benchmark - size: 1000, float SSE", 1000, 100);
SourceBench<tags::CPU::SSE, double > solver_bench_double_0_sse("Source   Benchmark - size: 1000, double SSE", 1000, 100);
#endif
