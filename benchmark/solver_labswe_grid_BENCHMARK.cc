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

#include <honei/lbm/solver_labswe_grid.hh>
#include <honei/swe/volume.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class LBMGSolverBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        LBMGSolverBench(const std::string & id, unsigned long size, int count) :
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

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLABSWEGrid<Tag_, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1., 1.5);

            solver.do_preprocessing();

            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(solver.solve());
            }
            LBMBenchmarkInfo benchinfo(SolverLABSWEGrid<tags::CPU, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>::get_benchmark_info(&grid, &data, &info));
            evaluate(benchinfo);
        }
};

LBMGSolverBench<tags::CPU, float> solver_bench_float_1("LBM Grid solver Benchmark - size: 250x250, float", 250, 25);
LBMGSolverBench<tags::CPU, double> solver_bench_double_1("LBM Grid solver Benchmark - size: 250x250, double", 250, 25);
LBMGSolverBench<tags::CPU::MultiCore, float> mc_solver_bench_float_1("MC LBM Grid solver Benchmark - size: 250x250, float", 250, 25);
LBMGSolverBench<tags::CPU::MultiCore, double> mc_solver_bench_double_1("MC LBM Grid solver Benchmark - size: 250x250, double", 250, 25);
#ifdef HONEI_SSE
LBMGSolverBench<tags::CPU::SSE, float> sse_solver_bench_float_1("SSE LBM Grid solver Benchmark - size: 250x250, float", 250, 25);
LBMGSolverBench<tags::CPU::SSE, double> sse_solver_bench_double_1("SSE LBM Grid solver Benchmark - size: 250x250, double", 250, 25);
LBMGSolverBench<tags::CPU::MultiCore::SSE, float> mcsse_solver_bench_float_1("MC SSE LBM Grid solver Benchmark - size: 250x250, float", 250, 25);
LBMGSolverBench<tags::CPU::MultiCore::SSE, double> mcsse_solver_bench_double_1("MC SSE LBM Grid solver Benchmark - size: 250x250, double", 250, 25);
#endif
#ifdef HONEI_CUDA
LBMGSolverBench<tags::GPU::CUDA, float> cuda_solver_bench_float_1("CUDA LBM Grid solver Benchmark - size: 250x250, float", 250, 25);
#endif
