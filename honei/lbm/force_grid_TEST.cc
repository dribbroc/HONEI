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
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/lbm/force_grid.hh>
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
class ForceSlopeLBMGridTest :
    public TaggedTest<Tag_>
{
    public:
        ForceSlopeLBMGridTest(const std::string & type) :
            TaggedTest<Tag_>("force_slope_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(10);
            unsigned long g_w(10);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(5, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);


            for(unsigned long i(0) ; i < 9 ; ++i)
            {
                (*data.distribution_x)[i] = 1;
                (*data.distribution_y)[i] = 1;
            }

            ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(info, data, DataType_(9.81), grid.d_x, grid.d_y, grid.d_t, DataType_(0.01));

            DenseVector<DataType_> temp_1_ref(data.f_temp_1->copy());
            DenseVector<DataType_> temp_2_ref(data.f_temp_2->copy());
            DenseVector<DataType_> temp_3_ref(data.f_temp_3->copy());
            DenseVector<DataType_> temp_4_ref(data.f_temp_4->copy());
            DenseVector<DataType_> temp_5_ref(data.f_temp_5->copy());
            DenseVector<DataType_> temp_6_ref(data.f_temp_6->copy());
            DenseVector<DataType_> temp_7_ref(data.f_temp_7->copy());
            DenseVector<DataType_> temp_8_ref(data.f_temp_8->copy());

            //----------------------------------------------------
            Grid<D2Q9, DataType_> grid_2;

            ScenarioCollection::get_scenario(5, g_h, g_w, grid_2);

            PackedGridData<D2Q9, DataType_>  data_2;
            PackedGridInfo<D2Q9> info_2;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_2, info_2, data_2);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info_2, data_2);

            data_2.distribution_x->lock(lm_write_only);
            data_2.distribution_y->lock(lm_write_only);
            for(unsigned long i(0) ; i < 9 ; ++i)
            {
                (*data_2.distribution_x)[i] = 1;
                (*data_2.distribution_y)[i] = 1;
            }
            data_2.distribution_x->unlock(lm_write_only);
            data_2.distribution_y->unlock(lm_write_only);

            ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(info_2, data_2, DataType_(9.81), grid_2.d_x, grid_2.d_y, grid_2.d_t, DataType_(0.01));

            data_2.f_temp_1->lock(lm_read_only);
            data_2.f_temp_2->lock(lm_read_only);
            data_2.f_temp_3->lock(lm_read_only);
            data_2.f_temp_4->lock(lm_read_only);
            data_2.f_temp_5->lock(lm_read_only);
            data_2.f_temp_6->lock(lm_read_only);
            data_2.f_temp_7->lock(lm_read_only);
            data_2.f_temp_8->lock(lm_read_only);

            //std::cout << temp_2_ref << std::endl;
            //std::cout << *data_2.f_temp_2 << std::endl;
            for(unsigned long i(0) ; i < temp_2_ref.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_1)[i], temp_1_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_2)[i], temp_2_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_3)[i], temp_3_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_4)[i], temp_4_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_5)[i], temp_5_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_6)[i], temp_6_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_7)[i], temp_7_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_8)[i], temp_8_ref[i], std::numeric_limits<DataType_>::epsilon());
            }
            data_2.f_temp_1->unlock(lm_read_only);
            data_2.f_temp_2->unlock(lm_read_only);
            data_2.f_temp_3->unlock(lm_read_only);
            data_2.f_temp_4->unlock(lm_read_only);
            data_2.f_temp_5->unlock(lm_read_only);
            data_2.f_temp_6->unlock(lm_read_only);
            data_2.f_temp_7->unlock(lm_read_only);
            data_2.f_temp_8->unlock(lm_read_only);

            /*data.f_temp_1->lock(lm_read_only);
            std::cout << *data.f_temp_1 << std::endl;
            data.f_temp_1->unlock(lm_read_only);
            data_2.f_temp_1->lock(lm_read_only);
            std::cout << *data_2.f_temp_1 << std::endl;
            data_2.f_temp_1->unlock(lm_read_only);*/
        }

};

#ifdef HONEI_SSE
ForceSlopeLBMGridTest<tags::CPU::SSE, float> force_slope_test_float_sse("float");
ForceSlopeLBMGridTest<tags::CPU::SSE, double> force_slope_test_double("double");
#endif

#ifdef HONEI_CUDA
ForceSlopeLBMGridTest<tags::GPU::CUDA, float> force_slope_test_float_cuda("CUDA float");
#endif

template <typename Tag_, typename DataType_>
class ForceFrictionLBMGridTest :
    public TaggedTest<Tag_>
{
    public:
        ForceFrictionLBMGridTest(const std::string & type) :
            TaggedTest<Tag_>("force_friction_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(10);
            unsigned long g_w(10);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(5, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);


            for(unsigned long i(0) ; i < 9 ; ++i)
            {
                (*data.distribution_x)[i] = 1;
                (*data.distribution_y)[i] = 1;
            }

            ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(info, data, DataType_(9.81), grid.d_x, grid.d_y, grid.d_t, DataType_(0.01));

            DenseVector<DataType_> temp_1_ref(data.f_temp_1->copy());
            DenseVector<DataType_> temp_2_ref(data.f_temp_2->copy());
            DenseVector<DataType_> temp_3_ref(data.f_temp_3->copy());
            DenseVector<DataType_> temp_4_ref(data.f_temp_4->copy());
            DenseVector<DataType_> temp_5_ref(data.f_temp_5->copy());
            DenseVector<DataType_> temp_6_ref(data.f_temp_6->copy());
            DenseVector<DataType_> temp_7_ref(data.f_temp_7->copy());
            DenseVector<DataType_> temp_8_ref(data.f_temp_8->copy());

            //----------------------------------------------------
            Grid<D2Q9, DataType_> grid_2;

            ScenarioCollection::get_scenario(5, g_h, g_w, grid_2);

            PackedGridData<D2Q9, DataType_>  data_2;
            PackedGridInfo<D2Q9> info_2;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_2, info_2, data_2);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info_2, data_2);

            data_2.distribution_x->lock(lm_write_only);
            data_2.distribution_y->lock(lm_write_only);
            for(unsigned long i(0) ; i < 9 ; ++i)
            {
                (*data_2.distribution_x)[i] = 1;
                (*data_2.distribution_y)[i] = 1;
            }
            data_2.distribution_x->unlock(lm_write_only);
            data_2.distribution_y->unlock(lm_write_only);

            ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(info_2, data_2, DataType_(9.81), grid_2.d_x, grid_2.d_y, grid_2.d_t, DataType_(0.01));

            data_2.f_temp_1->lock(lm_read_only);
            data_2.f_temp_2->lock(lm_read_only);
            data_2.f_temp_3->lock(lm_read_only);
            data_2.f_temp_4->lock(lm_read_only);
            data_2.f_temp_5->lock(lm_read_only);
            data_2.f_temp_6->lock(lm_read_only);
            data_2.f_temp_7->lock(lm_read_only);
            data_2.f_temp_8->lock(lm_read_only);
            for(unsigned long i(0) ; i < temp_1_ref.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_1)[i], temp_1_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_2)[i], temp_2_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_3)[i], temp_3_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_4)[i], temp_4_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_5)[i], temp_5_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_6)[i], temp_6_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_7)[i], temp_7_ref[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data_2.f_temp_8)[i], temp_8_ref[i], std::numeric_limits<DataType_>::epsilon());
            }
            data_2.f_temp_1->unlock(lm_read_only);
            data_2.f_temp_2->unlock(lm_read_only);
            data_2.f_temp_3->unlock(lm_read_only);
            data_2.f_temp_4->unlock(lm_read_only);
            data_2.f_temp_5->unlock(lm_read_only);
            data_2.f_temp_6->unlock(lm_read_only);
            data_2.f_temp_7->unlock(lm_read_only);
            data_2.f_temp_8->unlock(lm_read_only);

            /*data.f_temp_1->lock(lm_read_only);
            std::cout << *data.f_temp_1 << std::endl;
            data.f_temp_1->unlock(lm_read_only);
            data_2.f_temp_1->lock(lm_read_only);
            std::cout << *data_2.f_temp_1 << std::endl;
            data_2.f_temp_1->unlock(lm_read_only);*/
        }

};

#ifdef HONEI_SSE
ForceFrictionLBMGridTest<tags::CPU::SSE, float> force_friction_test_float_sse("float");
ForceFrictionLBMGridTest<tags::CPU::SSE, double> force_friction_test_double("double");
#endif

#ifdef HONEI_CUDA
ForceFrictionLBMGridTest<tags::GPU::CUDA, float> force_friction_test_float_cuda("CUDA float");
#endif
