/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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
#include <honei/lbm/solver_labswe_grid.hh>
#include <honei/swe/post_processing.hh>
#include <honei/swe/volume.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/grid_partitioner.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

template <typename Tag_, typename DataType_>
class SolverLABSWEGridTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEGridTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            //Other matrices needed by solver:
            /// \todo
            DenseVector<DataType_> s_x(data.h->size(), DataType_(0.));
            DenseVector<DataType_> s_y(data.h->size(), DataType_(0.));
            DenseVector<DataType_> b(data.h->size(), DataType_(0.));

            SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1., g_w, g_h, &b);

            solver.set_source(&s_x, &s_y);
            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            std::cout << *grid.h << std::endl;
#endif
            TEST_CHECK(true);
        }

};
SolverLABSWEGridTest<tags::CPU, double> solver_test_double("double");

template <typename Tag_, typename DataType_>
class SolverLABSWEGridPartitionerTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEGridPartitionerTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_partitioner_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            std::vector<PackedGridInfo<D2Q9> > info_list;
            std::vector<PackedGridData<D2Q9, DataType_> > data_list;
            GridPartitioner<D2Q9, DataType_>::decompose(2, info, data, info_list, data_list);

            //Other matrices needed by solver:
            /// \todo
            DenseVector<DataType_> s_x_0(data_list[0].h->size(), DataType_(0.));
            DenseVector<DataType_> s_y_0(data_list[0].h->size(), DataType_(0.));
            DenseVector<DataType_> b_0(data_list[0].h->size(), DataType_(0.));

            DenseVector<DataType_> s_x_1(data_list[1].h->size(), DataType_(0.));
            DenseVector<DataType_> s_y_1(data_list[1].h->size(), DataType_(0.));
            DenseVector<DataType_> b_1(data_list[1].h->size(), DataType_(0.));

            SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver_0(&data_list[0], &info_list[0], 1., 1., 1., g_w, g_h, &b_0);

            SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver_1(&data_list[1], &info_list[1], 1., 1., 1., g_w, g_h, &b_1);

            solver_0.set_source(&s_x_0, &s_y_0);
            solver_1.set_source(&s_x_1, &s_y_1);
            solver_0.do_preprocessing();
            solver_1.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver_0.solve();
                solver_1.solve();

                GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list);
#ifdef SOLVER_POSTPROCESSING
                GridPartitioner<D2Q9, DataType_>::compose(info, data, info_list, data_list);
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE

            GridPartitioner<D2Q9, DataType_>::compose(info, data, info_list, data_list);
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            std::cout << *grid.h << std::endl;
#endif
            TEST_CHECK(true);
        }
};
//SolverLABSWEGridPartitionerTest<tags::CPU, double> solver_partitioner_test_double("double");
