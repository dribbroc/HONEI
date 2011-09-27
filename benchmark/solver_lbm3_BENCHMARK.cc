/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <string>
#endif

#include <honei/woolb3/solver_lbm3.hh>
#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/swe/volume.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class LBM3SolverBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        LBM3SolverBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            unsigned long g_h(_size);
            unsigned long g_w(_size);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.04), 15, 15);
            b1.value();

            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();

            Grid3<DataType_, 9> grid3(obstacles, h, b, u, v);
            PackedGrid3<DataType_, 9> pgrid3(grid3);
            SolverLBM3<Tag_, DataType_, 9, lbm::lbm_source_schemes::BED_FULL> solver3(grid3, pgrid3, 1., 1., 1., 1.5);
            solver3.do_preprocessing();


            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 25 ; ++j)
                        {
                            solver3.solve();
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            evaluate();
        }
};
LBM3SolverBench<tags::CPU, float> solver_bench_float_1("LBM 3 solver Benchmark - size: 1500x1500, float", 250, 5);
LBM3SolverBench<tags::CPU, double> solver_bench_double_1("LBM 3 solver Benchmark - size: 1500x1500, double", 250, 5);
