/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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

#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/solver_lbm_fsi.hh>
#include <honei/lbm/scan_conversion_fsi.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/solid.hh>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/backends/cuda/operations.hh>

using namespace std;
using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class SolverLBMFSIFSIBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SolverLBMFSIFSIBench(const std::string & id, unsigned long size, int count) :
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


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedSolidData<D2Q9, DataType_>  solids;
            PackedGridInfo<D2Q9> info;

            DenseMatrix<bool> line(g_h, g_w, false);
            DenseMatrix<bool> bound(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            DenseMatrix<bool> sol(g_h, g_w, false);

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid, data, solids, line, bound, stf, sol, *grid.obstacles);

            SolverLBMFSI<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, &solids, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            //Directly dealing with omega-coordinates
            Line<DataType_, lbm_solid_dims::D2> line_1_0(DataType_(5) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(20) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_2_0(DataType_(10)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_3_0(DataType_(10)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_4_0(DataType_(5)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(20) * grid.d_y);

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(4);
            tri_0.add_line(line_1_0);
            tri_0.add_line(line_2_0);
            tri_0.add_line(line_3_0);
            tri_0.add_line(line_4_0);
            tri_0.value();

            ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_0, true);

            solids.current_u = DataType_(1./2.* grid.d_x);
            solids.current_v = DataType_(0.);
            solver.do_preprocessing();

            for(int i = 0; i < _count; ++i)
            {
                if (i < 80)
                {
                    Line<DataType_, lbm_solid_dims::D2> line_1_i(DataType_(5. + i/2.) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10. + i/2.)* grid.d_x, DataType_(20) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_2_i(DataType_(10. + i/2.)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10. + i/2.)* grid.d_x, DataType_(25) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_3_i(DataType_(10. + i/2.)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5. + i/2.)* grid.d_x, DataType_(25) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_4_i(DataType_(5. + i/2.)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5. + i/2.)* grid.d_x, DataType_(20) * grid.d_y);

                    Polygon<DataType_, lbm_solid_dims::D2> tri_i(4);
                    tri_i.add_line(line_1_i);
                    tri_i.add_line(line_2_i);
                    tri_i.add_line(line_3_i);
                    tri_i.add_line(line_4_i);
                    tri_i.value();
                    ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_i, true);
                }
                else
                {
                    solids.current_u = DataType_(0);
                }
                BENCHMARK(
                        solver.solve(1ul);
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda_thread_synchronize();
                        );
            }
            evaluate();
            data.destroy();
            info.destroy();
            solids.destroy();
        }
};

SolverLBMFSIFSIBench<tags::CPU, float> collide_stream_grid_bench_float("SolverLBMFSIFSIBenchmark - size: 129, float", 129, 100);
SolverLBMFSIFSIBench<tags::CPU, double> collide_stream_grid_bench_double("SolverLBMFSIFSIBenchmark - size: 129, double", 129, 100);
#ifdef HONEI_CUDA
SolverLBMFSIFSIBench<tags::GPU::CUDA, float> collide_stream_grid_bench_float_cuda("SolverLBMFSIFSIBenchmark CUDA - size: 450*2, float", 450*2, 100);
#endif
