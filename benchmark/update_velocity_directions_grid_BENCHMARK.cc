/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
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

#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/backends/cuda/operations.hh>

using namespace std;
using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class UpdateVelocityDirectionsGridBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        UpdateVelocityDirectionsGridBench(const std::string & id, unsigned long size, int count) :
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
            grid.obstacles = new DenseMatrix<bool>(obstacles);
            grid.h = new DenseMatrix<DataType_>(h);
            grid.u = new DenseMatrix<DataType_>(u);
            grid.v = new DenseMatrix<DataType_>(v);
            grid.b = new DenseMatrix<DataType_>(b);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, 1., 1., 1., 1.5);

            solver.do_preprocessing();

            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                        (UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::
                         value(info, data));
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            evaluate();
            data.destroy();
            info.destroy();
            grid.destroy();
        }
};

UpdateVelocityDirectionsGridBench<tags::CPU, float> update_velocity_directions_grid_bench_float("UpdateVelocityDirectionsGridBench - size: 2000, float", 2000, 10);
UpdateVelocityDirectionsGridBench<tags::CPU, double> update_velocity_directions_grid_bench_double("UpdateVelocityDirectionsGridBench - size: 2000, double", 2000, 10);
#ifdef HONEI_SSE
UpdateVelocityDirectionsGridBench<tags::CPU::SSE, float> sse_update_velocity_directions_grid_bench_float("SSE UpdateVelocityDirectionsGridBench - size: 2000, float", 2000, 10);
UpdateVelocityDirectionsGridBench<tags::CPU::SSE, double> sse_colliupdate_velocity_directions_bench_double("SSE UpdateVelocityDirectionsGridBench - size: 2000, double", 2000, 10);
#endif
#ifdef HONEI_CUDA
UpdateVelocityDirectionsGridBench<tags::GPU::CUDA, float> cuda_update_velocity_directions_grid_bench_float("CUDA UpdateVelocityDirectionsGridBench - size: 2000, float", 2000, 10);
#endif
#ifdef HONEI_CELL
UpdateVelocityDirectionsGridBench<tags::Cell, float> cell_update_velocity_directions_grid_bench_float("Cell UpdateVelocityDirectionsGridBench - size: 250, float", 250, 25);
#endif
