/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@mathematik.uni-dortmund.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/util/kpnetcdffile.hh>
#include <honei/la/difference.hh>

#ifdef DEBUG
#define SOLVER_VERBOSE
#endif

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLBMGridTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;
            grid.d_x = 0.01;
            grid.d_y = 0.01;
            grid.d_t = 0.01;
            grid.tau = 1.1;

            std::string path(HONEI_SOURCEDIR);
            path += "/honei/lbm/sample.nc";
            KPInitialConditions cond;
            std::tr1::shared_ptr<Field> b;
            std::tr1::shared_ptr<Field> u1;
            std::tr1::shared_ptr<Field> u2;
            std::tr1::shared_ptr<Field> u3;
            KPNetCDFFile file(path, cond, b, u1, u2, u3);
            unsigned long g_h = cond.getNx();
            unsigned long g_w = cond.getNy();
            grid.obstacles = new DenseMatrix<bool>(g_h, g_w, false);
            grid.b = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0));
            grid.h = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0));
            grid.u = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0));
            grid.v = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0));
            for (unsigned long i(0) ; i < g_h * g_w ; ++i)
            {
                (*grid.b)(i/g_w, i%g_w) = b->data[i];
                (*grid.h)(i/g_w, i%g_w) = u1->data[i];
                (*grid.u)(i/g_w, i%g_w) = u2->data[i];
                (*grid.v)(i/g_w, i%g_w) = u3->data[i];
            }

            Difference<tags::CPU>::value(*grid.h, *grid.b);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif
            for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));
        }

};

//SolverLBMGridTest<tags::CPU, float> solver_test_float("float");
/*SolverLBMGridTest<tags::CPU, double> solver_test_double("double");
SolverLBMGridTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLBMGridTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");
#ifdef HONEI_SSE
SolverLBMGridTest<tags::CPU::SSE, float> sse_solver_test_float("float");
SolverLBMGridTest<tags::CPU::SSE, double> sse_solver_test_double("double");
SolverLBMGridTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_test_float("float");
SolverLBMGridTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLBMGridTest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
SolverLBMGridTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_solver_test_float("float");
#endif
#ifdef HONEI_CELL
SolverLBMGridTest<tags::Cell, float> cell_solver_test_float("float");
#endif
*/
