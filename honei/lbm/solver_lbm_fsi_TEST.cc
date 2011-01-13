/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
#include <honei/lbm/solver_lbm_fsi.hh>
#include <honei/lbm/scan_conversion_fsi.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLBMFSITest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSITest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(0, g_h, g_w, grid);

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

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
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

            info.destroy();
            data.destroy();
            grid.destroy();
            solids.destroy();
        }

};
SolverLBMFSITest<tags::CPU, float> solver_test_float("float");
#ifdef HONEI_CUDA
SolverLBMFSITest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
#endif


template <typename Tag_, typename DataType_>
class SolverLBMFSIStationaryTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSIStationaryTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_stationary_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(0, g_h, g_w, grid);

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

            ScanConversionFSI<Tag_>::value(grid, info, data, solids, tri_0, true);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
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
            info.destroy();
            data.destroy();
            grid.destroy();
            solids.destroy();
        }

};
//SolverLBMFSIStationaryTest<tags::CPU, float> solver_test_stat_float("float");

template <typename Tag_, typename DataType_>
class SolverLBMFSINonStationaryTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSINonStationaryTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_stationary_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(100);
            unsigned long g_w(100);
            unsigned long timesteps(100);


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
            solids.current_v = -DataType_(1./2.* grid.d_x);//DataType_(0.);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                //Directly dealing with omega-coordinates
                if (i < 80)
                {
                    Line<DataType_, lbm_solid_dims::D2> line_1_i(DataType_(5.+ i/2.) * grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_2_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_3_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_4_i(DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);

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
                solver.solve(4ul);
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
        }

};
//SolverLBMFSINonStationaryTest<tags::CPU, float> solver_test_nstat_float("float");
#ifdef HONEI_CUDA
//SolverLBMFSINonStationaryTest<tags::GPU::CUDA, float> solver_test_nstat_cuda_float(" CUDA float");
#endif

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_1 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_1(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_1<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 50);
            unsigned long g_w(_discr_lev * 50);
            unsigned long timesteps(100000);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(9, g_h, g_w, grid);
            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

            Cuboid<DataType_> a(*grid.h, 5*_discr_lev, 5*_discr_lev, DataType_(0.02), g_w/2, g_h/2);
            a.value();

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

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(5.* DataType_(_discr_lev) * grid.d_x * 5. * DataType_(_discr_lev) * grid.d_y * 0.02 + ((g_w*grid.d_x) * (g_h*grid.d_x) * 0.05));
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            grid.h->lock(lm_read_only);
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.h->unlock(lm_read_only);
            std::cout << "Vol.: " << vol << std::endl;
        }

};
//SolverLBMFSITest_MC_1<tags::GPU::CUDA, float> solver_test_float_mc1("float", 1ul);

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_2 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_2(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_2<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 50);
            unsigned long g_w(_discr_lev * 50);
            unsigned long timesteps(100000);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(9, g_h, g_w, grid);
            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

            Cuboid<bool> d1(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(5*_discr_lev), true, (signed long)(g_w/2), 0l);
            d1.value();
            Cuboid<bool> d2(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(5*_discr_lev), true, (signed long)g_w/2, (signed long)(g_h - (20 * _discr_lev)));
            d2.value();
            Cuboid<DataType_> a(*grid.h, 50*_discr_lev, 20*_discr_lev, DataType_(0.02), g_w - (20 * _discr_lev), 0l);
            a.value();

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

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 100, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(20.* DataType_(_discr_lev) * grid.d_x * 50. * DataType_(_discr_lev) * grid.d_y * 0.02 + ((g_w*grid.d_x) * (g_h*grid.d_x) * 0.05) - 2.*(5. * DataType_(_discr_lev) * grid.d_x * 20. * DataType_(_discr_lev) * grid.d_y* 0.05));
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            grid.h->lock(lm_read_only);
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.h->unlock(lm_read_only);
            std::cout << "Vol.: " << vol << std::endl;
        }

};
//SolverLBMFSITest_MC_2<tags::GPU::CUDA, float> solver_test_float_mc2("float", 1ul);

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_3 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_3(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_3<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 50);
            unsigned long g_w(_discr_lev * 50);
            unsigned long timesteps(100000);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(9, g_h, g_w, grid);
            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

            Cuboid<DataType_> a(*grid.h, 5*_discr_lev, 5*_discr_lev, DataType_(0.02), g_w/2, g_h/2);
            a.value();

            //build up the hill:
            for(unsigned long i(0) ; i <  g_h ; ++i)
            {
                for(unsigned long j(0) ; j <  g_w ; ++j)
                {
                    double x(j *  grid.d_x);
                    double y(i *  grid.d_y);
                    //if(sqrt(y * y + x * x) >= 0.4)
                    (*grid.b)(i , j) = 0.0000001 * (x*x +y*y);
                }
            }

            for(unsigned long i(0) ; i <  g_h ; ++i)
            {
                for(unsigned long j(0) ; j <  g_w ; ++j)
                {
                    (*grid.h)(i , j) -= (*grid.b)(i , j);
                }
            }
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

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 100, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(5.* DataType_(_discr_lev) * grid.d_x * 5. * DataType_(_discr_lev) * grid.d_y * 0.02 + ((g_w*grid.d_x) * (g_h*grid.d_x) * 0.05));
            DataType_ ana_vol_hill(10e-7 * 2. * (50 * 50 * 50 * 50) / 3.);
            std::cout << "Analytical Vol.: " << ana_vol - ana_vol_hill << " ";
            grid.h->lock(lm_read_only);
            grid.b->lock(lm_read_only);
            DataType_ vol_h = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            DataType_ vol_b = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.b, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.b->unlock(lm_read_only);
            std::cout << "Vol.: " << vol_h - vol_b << std::endl;
        }

};
//SolverLBMFSITest_MC_3<tags::GPU::CUDA, float> solver_test_float_mc3("float", 1ul);

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_4 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_4(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_4<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 50);
            unsigned long g_w(_discr_lev * 50);
            unsigned long timesteps(100000);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(9, g_h, g_w, grid);
            grid.d_x = 0.01;
            grid.d_y = 0.01;
            grid.d_t = 0.01;
            grid.tau = 1.1;
            for(unsigned long i(0) ; i < g_h ; ++i)
                for(unsigned long j(0) ; j < g_w ; ++j)
                    (*grid.h)[i][j] = DataType_(0);

            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

            Cuboid<bool> d1(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(5*_discr_lev), true, (signed long)(g_w/2), 0l);
            d1.value();
            Cuboid<bool> d2(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(5*_discr_lev), true, (signed long)g_w/2, (signed long)(g_h - (20 * _discr_lev)));
            d2.value();
            Cuboid<DataType_> a(*grid.h, 50*_discr_lev, 20*_discr_lev, DataType_(0.02), g_w - (20 * _discr_lev), 0l);
            a.value();

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

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 100, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(20.* DataType_(_discr_lev) * grid.d_x * 50. * DataType_(_discr_lev) * grid.d_y * 0.02);
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            grid.h->lock(lm_read_only);
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.h->unlock(lm_read_only);
            std::cout << "Vol.: " << vol << std::endl;
        }

};
//SolverLBMFSITest_MC_4<tags::GPU::CUDA, float> solver_test_float_mc4("float", 1ul);

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_5 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_5(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_5<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 100);
            unsigned long g_w(_discr_lev * 100);
            unsigned long timesteps(100000);

            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);
            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

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
            solids.current_v = -DataType_(1./2.* grid.d_x);//DataType_(0.);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                if (i < 30)
                {
                    Line<DataType_, lbm_solid_dims::D2> line_1_i(DataType_(5.+ i/2.) * grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_2_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_3_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_4_i(DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);

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
                    solids.current_v = DataType_(0);
                }
                solver.solve(2ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(((g_w*grid.d_x) * (g_h*grid.d_x) * 0.05) - (10. * DataType_(_discr_lev) * grid.d_x * 10. * DataType_(_discr_lev) * grid.d_y * 0.05));
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            grid.h->lock(lm_read_only);
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.h->unlock(lm_read_only);
            std::cout << "Vol.: " << vol << std::endl;
        }

};
//SolverLBMFSITest_MC_5<tags::GPU::CUDA, float> solver_test_float_mc5("float", 1ul);

template <typename Tag_, typename DataType_>
class SolverLBMFSITest_MC_6 :
    public TaggedTest<Tag_>
{
    private:
        unsigned long _discr_lev;
    public:
        SolverLBMFSITest_MC_6(const std::string & type, unsigned long discr_lev) :
            TaggedTest<Tag_>("solver_lbm_fsi_test_MC_6<" + type + ">")
        {
            _discr_lev = discr_lev;
        }

        virtual void run() const
        {
            unsigned long g_h(_discr_lev * 100);
            unsigned long g_w(_discr_lev * 100);
            unsigned long timesteps(100000);

            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);
            grid.d_x /= DataType_(_discr_lev);
            grid.d_y /= DataType_(_discr_lev);
            grid.d_t /= DataType_(_discr_lev);

            Cuboid<DataType_> a(*grid.h, 15*_discr_lev, 15*_discr_lev, DataType_(0.02), g_w/2, g_h/2);
            a.value();

            /*Cuboid<bool> b(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(20*_discr_lev), true, 0, (signed long)(g_h/2));
            b.value();
            Cuboid<bool> b1(*grid.obstacles, (signed long)(20*_discr_lev), (signed long)(20*_discr_lev), true, (signed long)(g_h/2), 0);
            b1.value();*/
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
            solids.current_v = -DataType_(1./2.* grid.d_x);//DataType_(0.);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                if (i < 80)
                {
                    Line<DataType_, lbm_solid_dims::D2> line_1_i(DataType_(5.+ i/2.) * grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_2_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_3_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_4_i(DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);

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
                    solids.current_v = -DataType_(0);
                }
                solver.solve(4ul);
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
                Difference<Tag_>::value(last_u, *data.u);
                last_u.lock(lm_read_only);
                DataType_ n(Norm<vnt_l_two, false, tags::CPU>::value(last_u));
                last_u.unlock(lm_read_only);
                if( n <= std::numeric_limits<DataType_>::epsilon())
                    break;

                copy<Tag_>(*data.u, last_u);
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif

            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(15.* DataType_(_discr_lev) * grid.d_x * 15. * DataType_(_discr_lev) * grid.d_y * 0.02 + ((g_w*grid.d_x) * (g_h*grid.d_x) * 0.05) - (10. * DataType_(_discr_lev) * grid.d_x * 10. * DataType_(_discr_lev) * grid.d_y * 0.05) - 0.*(20. * DataType_(_discr_lev) * grid.d_x * 20. * DataType_(_discr_lev) * grid.d_y * 0.05));
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            grid.h->lock(lm_read_only);
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(*grid.h, DataType_(0), DataType_(g_w*grid.d_x), grid.d_x, grid.d_y);
            grid.h->unlock(lm_read_only);
            std::cout << "Vol.: " << vol << std::endl;
        }

};
//SolverLBMFSITest_MC_6<tags::GPU::CUDA, float> solver_test_float_mc6("float", 1ul);

