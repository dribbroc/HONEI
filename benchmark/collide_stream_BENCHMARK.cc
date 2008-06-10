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

#include <honei/lbm/collide_stream.hh>

using namespace std;
using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_, typename Dir_>
class CollideStreamBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        CollideStreamBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseMatrix<DataType_> test(_size, _size, DataType_(1.23456));
            DataType_ e_x(1.);
            DataType_ e_y(1.);
            DataType_ tau(1.);

            DenseMatrix<DataType_> result(_size, _size);

            for(int i = 0; i < _count; ++i)
            {

                BENCHMARK((CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, Dir_>::value(result, test, test, test, test, e_x, e_y, tau)));

            }
            evaluate();
        }
};

CollideStreamBench<tags::CPU, float, D2Q9::DIR_0> solver_bench_float_0("CollideStream DIR_0 Benchmark - size: 1000, float", 1000, 100);
CollideStreamBench<tags::CPU, double, D2Q9::DIR_0> solver_bench_double_0("CollideStream DIR_0 Benchmark - size: 1000, double", 1000, 100);
CollideStreamBench<tags::CPU, float, D2Q9::DIR_1> solver_bench_float_1("CollideStream DIR_1 Benchmark - size: 1000, float", 1000, 100);
CollideStreamBench<tags::CPU, double, D2Q9::DIR_1> solver_bench_double_1("CollideStream DIR_1 Benchmark - size: 1000, double", 1000, 100);
CollideStreamBench<tags::CPU, float, D2Q9::DIR_2> solver_bench_float_2("CollideStream DIR_0 Benchmark - size: 1000, float", 1000, 100);
CollideStreamBench<tags::CPU, double, D2Q9::DIR_2> solver_bench_double_2("CollideStream DIR_0 Benchmark - size: 1000, double", 1000, 100);
