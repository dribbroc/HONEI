/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. HONEI is free software;
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

#include <honei/lbm/solver_lbm_fsi.hh>
#include <honei/lbm/scan_conversion_fsi.hh>
//#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/lbm/bitmap_io.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/la/difference.hh>
#include <honei/la/sum.hh>
#include <honei/la/scale.hh>
#include <honei/math/matrix_io.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

#define SOLVER_VERBOSE 1
#define SOLVER_POSTPROCESSING 1


int main(int argc, char ** argv)
{
    float dx(15);
    unsigned long _max_timesteps(5000);
    if (argc != 3)
    {
        std::cout<<"Usage: norwegen delta_t tau"<<std::endl;
        exit(1);
    }
    float _dt(atof(argv[1]));
    float _tau(atof(argv[2]));

    SparseMatrixELL<float> pre_h_b(MatrixIO<io_formats::ELL>::read_matrix("malp_initial.ell", float(0)));
    DenseMatrix<float> h_b(pre_h_b.rows(), pre_h_b.columns());

    SparseMatrixELL<float> pre_pre_s(MatrixIO<io_formats::ELL>::read_matrix("malp_domain.ell", float(0)));
    DenseMatrix<float> pre_s(pre_h_b.rows(), pre_h_b.columns());

    SparseMatrixELL<float> pre_b(MatrixIO<io_formats::ELL>::read_matrix("malp_bottom.ell", float(0)));
    DenseMatrix<float> b(pre_h_b.rows(), pre_h_b.columns());

    unsigned long height(pre_h_b.rows());
    unsigned long width(pre_h_b.columns());
    DenseMatrix<float> h((unsigned long)height, (unsigned long)width, float(0));
    DenseMatrix<float> u((unsigned long)height, (unsigned long)width, float(0));
    DenseMatrix<float> v((unsigned long)height, (unsigned long)width, float(0));
    DenseMatrix<bool> o((unsigned long)height, (unsigned long)width, float(0));
    DenseMatrix<bool> s((unsigned long)height, (unsigned long)width, float(0));

    for (unsigned long i = 0; i < height; ++i)
    {
        for (unsigned long j = 0; j < width; ++j)
        {
            h_b[i][j] = pre_h_b(i , j);
            pre_s[i][j] = pre_pre_s(i , j);
            s[i][j] = pre_s[i][j] == 110 ? true : false;
            b[i][j] = s(i , j) ? float(0) : pre_s(i , j);
            h[i][j] = h_b[i][j] - b[i][j];
        }
    }
    //PostProcessing<GNUPLOT>::value(b, 1, b.columns(), b.rows(), 1000);
    //PostProcessing<GNUPLOT>::value(h, 1, b.columns(), b.rows(), 2000);

    DenseMatrix<float> twenty((unsigned long)height, (unsigned long)width, float(20));
    Sum<tags::CPU>::value(b ,twenty);
    Sum<tags::CPU>::value(h ,twenty);

    Scale<tags::CPU>::value(b, float(1./1000));
    Scale<tags::CPU>::value(h, float(1./1000));
    dx /= 1000;

    Grid<D2Q9, float> grid;
    grid.h = &h;
    grid.u = &u;
    grid.v = &v;
    grid.b = &b;
    grid.obstacles = &o;

    grid.d_x = dx;
    grid.d_y = dx;
    grid.d_t = _dt;
    grid.tau = _tau;

    PackedGridData<D2Q9, float>  data;
    PackedSolidData<D2Q9, float>  solids;
    PackedGridInfo<D2Q9> info;

    GridPacker<D2Q9, NOSLIP, float>::pack(grid, info, data);
    GridPacker<D2Q9, lbm_boundary_types::NOSLIP, float>::cuda_pack(info, data);
    GridPackerFSI<D2Q9, NOSLIP, float>::allocate(data, solids);

    DenseMatrix<bool> line(s.copy());
    DenseMatrix<bool> bound(s.copy());
    DenseMatrix<bool> stf(s.copy());
    GridPackerFSI<D2Q9, NOSLIP, float>::pack(grid, data, solids, line, bound, stf, s, *grid.obstacles);

    SolverLBMFSI<tags::CPU, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, &solids, grid.d_x, grid.d_y, grid.d_t, grid.tau);

    solids.current_u = float(0);
    solids.current_v = float(0);
    solver.do_preprocessing();

    for(unsigned long i(0); i < _max_timesteps; ++i)
    {
#ifdef SOLVER_VERBOSE
        std::cout<<"Timestep: " << i << "/" << _max_timesteps << std::endl;
#endif
        solver.solve(0ul);
#ifdef SOLVER_POSTPROCESSING
        //PostProcessing<GNUPLOT>::value(*grid.h, _max_timesteps - 1 , h.columns(), h.rows(), i);
        //Sum<tags::CPU>::value(*grid.h, *grid.b);
#endif
    }
    solver.do_postprocessing();
    GridPacker<D2Q9, NOSLIP, float>::unpack(grid, info, data);
    std::string filename("out_");
    filename += stringify(_dt);
    filename += "_";
    filename += stringify(_tau);
    filename +=".dat";
    PostProcessing<GNUPLOT>::value(*grid.h, 1, h.columns(), h.rows(), 0, filename);
    //GridPacker<D2Q9, NOSLIP, float>::unpack(grid, info, data);

    DenseMatrix<float> result(grid.h->copy());
    Sum<tags::CPU>::value(result, b);
    //std::cout << result << std::endl;
}
