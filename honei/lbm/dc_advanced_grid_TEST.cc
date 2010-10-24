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
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/dc_util.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;
using namespace lbm;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING
//#define DRIVEN_CAVITY_OUTPUT_TESTLINE
//#define DRIVEN_CAVITY_OUTPUT_ACCURACY

template <typename Tag_, typename DataType_>
class SolverLABNAVSTOGridDCTest :
    public TaggedTest<Tag_>
{
    private:
        DataType_ _dx, _dy, _dt, _reynolds, _tau, _U;
    public:
        SolverLABNAVSTOGridDCTest(const std::string & type, DataType_ dx, DataType_ dy, DataType_ dt, DataType_ reynolds, DataType_ tau, DataType_ U) :
            TaggedTest<Tag_>("solver_lbm_grid_dc_advanced_test<" + type + ">")
        {
            _dx = dx;
            _dy = dy;
            _dt = dt;
            _reynolds = reynolds;
            _tau = tau;
            _U = U;
        }

        virtual void run() const
        {
            unsigned long g_h(129);
            unsigned long g_w(129);
            unsigned long timesteps(100000);
            DataType_ dx(_dx);
            DataType_ dy(_dy);
            DataType_ dt(_dt);
            DataType_ tau(_tau);
            DataType_ lid_U(Reynolds::adjust_veloc(double(_reynolds), double(0.), g_w, dx, dt, tau));

            std::cout << "U: " << lid_U << std::endl;
            std::cout << "Reynolds: " << Reynolds::value(lid_U, g_w, dx, dt, tau) << std::endl;

            Grid<D2Q9, DataType_> grid;

            //just take some scenario and adjust it;
            ScenarioCollection::get_scenario(8, g_h, g_w, grid);
            grid.d_x = dx;
            grid.d_y = dy;
            grid.d_t = dt;
            grid.tau = tau;

            //init lid velocity
            for (unsigned long i(0) ; i < g_w ; ++i)
            {
                (*grid.u)[0][i] = lid_U;
            }

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);


            SolverLBMGrid<Tag_, lbm_applications::LABNAVSTO, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            std::cout << "dx " << grid.d_x << std::endl;
            std::cout << "dy " << grid.d_y << std::endl;
            std::cout << "dt " << grid.d_t << std::endl;
            std::cout << "tau " << grid.tau << std::endl;

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            DenseVector<DataType_> last_u(data.u->copy());
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                //after solving, reset lid velocity by using h_index
                for (unsigned long j(0) ; j < g_w ; ++j)
                {
                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = lid_U;
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = DataType_(0.);

                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);

                    if(j > 0)
                    {
                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);

                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                    }

                }

                //solver.solve();
                EquilibriumDistributionGrid<Tag_, lbm_applications::LABNAVSTO>::
                    value(DataType_(9.81), (grid.d_x / grid.d_t) * (grid.d_x / grid.d_t), info, data);

                CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                    value(info,
                          data,
                          tau);
                ///Boundary correction:
                UpdateVelocityDirectionsGrid<Tag_, lbm_boundary_types::NOSLIP>::
                    value(info, data);

                //extract velocities out of h from previous timestep:
                ExtractionGrid<Tag_, lbm_modes::WET>::value(info, data, DataType_(10e-5));

                //after solving, reset lid velocity by using h_index
                for (unsigned long j(0) ; j < g_w ; ++j)
                {
                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = lid_U;
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = DataType_(0.);

                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);

                    if(j > 0)
                    {
                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);

                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                    }

                }
                std::cout << i << std::endl;
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack_u(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.u, 999, g_w, g_h, i);
#endif
                if(Norm<vnt_l_two, false, Tag_>::value(Difference<Tag_>::value(last_u, *data.u)) <= std::numeric_limits<DataType_>::epsilon())
                    break;

                last_u = data.u->copy();
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack_u(grid, info, data);
            DenseVector<DataType_> test_line(g_h);
            std::cout<<"Index: " << (g_w - 1) / 2 << std::endl;
            for(unsigned long i(0); i < g_h; ++i)
            {
                test_line[i] = DataType_( (*grid.u)(i,( g_w - 1 ) / 2)/ lid_U);
            }
#ifdef DRIVEN_CAVITY_OUTPUT_TESTLINE
            std::string filename;
            std::ofstream file;
            filename = "out_lbm_navsto_dc_relative_grid.dat";
            file.open(filename.c_str());
            for(unsigned long i(0); i < g_h; ++i)
            {
                file << stringify(i) + " " + stringify(test_line[i]) + "\n";
            }
            file.close();
#endif
            std::cout<<"Result:"<<test_line<<std::endl;

            //Reference data by Ghia et al. 1982 for reynolds number of 100:
            DenseVector<double> ref_result_100(17);
            unsigned long indices_100[17];

            ref_result_100[0] = double(1);
            indices_100[0] = 0;

            ref_result_100[1] = double(0.84123);
            indices_100[1] = 3;

            ref_result_100[2] = double(0.78871);
            indices_100[2] = 4;

            ref_result_100[3] = double(0.73722);
            indices_100[3] = 5;

            ref_result_100[4] = double(0.68717);
            indices_100[4] = 6;

            ref_result_100[5] = double(0.23151);
            indices_100[5] = 19;

            ref_result_100[6] = double(0.00332);
            indices_100[6] = 34;

            ref_result_100[7] = double(-0.13641);
            indices_100[7] = 49;

            ref_result_100[8] = double(-0.20581);
            indices_100[8] = 64;

            ref_result_100[9] = double(-0.21090);
            indices_100[9] = 70;

            ref_result_100[10] = double(-0.15662);
            indices_100[10] = 92;

            ref_result_100[11] = double(-0.10150);
            indices_100[11] = 106;

            ref_result_100[12] = double(-0.06434);
            indices_100[12] = 115;

            ref_result_100[13] = double(-0.04775);
            indices_100[13] = 119;

            ref_result_100[14] = double(-0.04192);
            indices_100[14] = 120;

            ref_result_100[15] = double(-0.03717);
            indices_100[15] = 121;

            ref_result_100[16] = double(0.);
            indices_100[16] = 128;

            DenseVector<DataType_> diff(test_line.copy());

            for(unsigned long i(0); i < 17; ++i)
            {
                diff[indices_100[i]] = ref_result_100[i];
                //TEST_CHECK_EQUAL_WITHIN_EPS(test_line[indices_100[i]], ref_result_100[i], 0.5);
            }

            //compare with HUEBNER
            DenseVector<DataType_> ref_result_100_hubner(g_w);
            unsigned long i(0);

            ref_result_100_hubner[i] = 1.000000e+00;
            ++i;
            ref_result_100_hubner[i] = 9.482407e-01;
            ++i;
            ref_result_100_hubner[i] = 8.958706e-01;
            ++i;
            ref_result_100_hubner[i] = 8.434837e-01;
            ++i;
            ref_result_100_hubner[i] = 7.916057e-01;
            ++i;
            ref_result_100_hubner[i] = 7.407020e-01;
            ++i;
            ref_result_100_hubner[i] = 6.911695e-01;
            ++i;
            ref_result_100_hubner[i] = 6.433326e-01;
            ++i;
            ref_result_100_hubner[i] = 5.974428e-01;
            ++i;
            ref_result_100_hubner[i] = 5.536819e-01;
            ++i;
            ref_result_100_hubner[i] = 5.121672e-01;
            ++i;
            ref_result_100_hubner[i] = 4.729579e-01;
            ++i;
            ref_result_100_hubner[i] = 4.360633e-01;
            ++i;
            ref_result_100_hubner[i] = 4.014499e-01;
            ++i;
            ref_result_100_hubner[i] = 3.690501e-01;
            ++i;
            ref_result_100_hubner[i] = 3.387690e-01;
            ++i;
            ref_result_100_hubner[i] = 3.104923e-01;
            ++i;
            ref_result_100_hubner[i] = 2.840917e-01;
            ++i;
            ref_result_100_hubner[i] = 2.594308e-01;
            ++i;
            ref_result_100_hubner[i] = 2.363697e-01;
            ++i;
            ref_result_100_hubner[i] = 2.147690e-01;
            ++i;
            ref_result_100_hubner[i] = 1.944924e-01;
            ++i;
            ref_result_100_hubner[i] = 1.754096e-01;
            ++i;
            ref_result_100_hubner[i] = 1.573977e-01;
            ++i;
            ref_result_100_hubner[i] = 1.403425e-01;
            ++i;
            ref_result_100_hubner[i] = 1.241392e-01;
            ++i;
            ref_result_100_hubner[i] = 1.086928e-01;
            ++i;
            ref_result_100_hubner[i] = 9.391818e-02;
            ++i;
            ref_result_100_hubner[i] = 7.973991e-02;
            ++i;
            ref_result_100_hubner[i] = 6.609181e-02;
            ++i;
            ref_result_100_hubner[i] = 5.291659e-02;
            ++i;
            ref_result_100_hubner[i] = 4.016524e-02;
            ++i;
            ref_result_100_hubner[i] = 2.779642e-02;
            ++i;
            ref_result_100_hubner[i] = 1.577586e-02;
            ++i;
            ref_result_100_hubner[i] = 4.075649e-03;
            ++i;
            ref_result_100_hubner[i] = -7.326301e-03;
            ++i;
            ref_result_100_hubner[i] = -1.844692e-02;
            ++i;
            ref_result_100_hubner[i] = -2.929853e-02;
            ++i;
            ref_result_100_hubner[i] = -3.988940e-02;
            ++i;
            ref_result_100_hubner[i] = -5.022425e-02;
            ++i;
            ref_result_100_hubner[i] = -6.030470e-02;
            ++i;
            ref_result_100_hubner[i] = -7.012976e-02;
            ++i;
            ref_result_100_hubner[i] = -7.969618e-02;
            ++i;
            ref_result_100_hubner[i] = -8.899887e-02;
            ++i;
            ref_result_100_hubner[i] = -9.803122e-02;
            ++i;
            ref_result_100_hubner[i] = -1.067854e-01;
            ++i;
            ref_result_100_hubner[i] = -1.152528e-01;
            ++i;
            ref_result_100_hubner[i] = -1.234240e-01;
            ++i;
            ref_result_100_hubner[i] = -1.312894e-01;
            ++i;
            ref_result_100_hubner[i] = -1.388392e-01;
            ++i;
            ref_result_100_hubner[i] = -1.460636e-01;
            ++i;
            ref_result_100_hubner[i] = -1.529531e-01;
            ++i;
            ref_result_100_hubner[i] = -1.594988e-01;
            ++i;
            ref_result_100_hubner[i] = -1.656923e-01;
            ++i;
            ref_result_100_hubner[i] = -1.715258e-01;
            ++i;
            ref_result_100_hubner[i] = -1.769927e-01;
            ++i;
            ref_result_100_hubner[i] = -1.820873e-01;
            ++i;
            ref_result_100_hubner[i] = -1.868048e-01;
            ++i;
            ref_result_100_hubner[i] = -1.911417e-01;
            ++i;
            ref_result_100_hubner[i] = -1.950958e-01;
            ++i;
            ref_result_100_hubner[i] = -1.986661e-01;
            ++i;
            ref_result_100_hubner[i] = -2.018528e-01;
            ++i;
            ref_result_100_hubner[i] = -2.046576e-01;
            ++i;
            ref_result_100_hubner[i] = -2.070833e-01;
            ++i;
            ref_result_100_hubner[i] = -2.091342e-01;
            ++i;
            ref_result_100_hubner[i] = -2.108155e-01;
            ++i;
            ref_result_100_hubner[i] = -2.121341e-01;
            ++i;
            ref_result_100_hubner[i] = -2.130976e-01;
            ++i;
            ref_result_100_hubner[i] = -2.137150e-01;
            ++i;
            ref_result_100_hubner[i] = -2.139960e-01;
            ++i;
            ref_result_100_hubner[i] = -2.139515e-01;
            ++i;
            ref_result_100_hubner[i] = -2.135930e-01;
            ++i;
            ref_result_100_hubner[i] = -2.129330e-01;
            ++i;
            ref_result_100_hubner[i] = -2.119843e-01;
            ++i;
            ref_result_100_hubner[i] = -2.107604e-01;
            ++i;
            ref_result_100_hubner[i] = -2.092753e-01;
            ++i;
            ref_result_100_hubner[i] = -2.075432e-01;
            ++i;
            ref_result_100_hubner[i] = -2.055785e-01;
            ++i;
            ref_result_100_hubner[i] = -2.033957e-01;
            ++i;
            ref_result_100_hubner[i] = -2.010094e-01;
            ++i;
            ref_result_100_hubner[i] = -1.984341e-01;
            ++i;
            ref_result_100_hubner[i] = -1.956842e-01;
            ++i;
            ref_result_100_hubner[i] = -1.927736e-01;
            ++i;
            ref_result_100_hubner[i] = -1.897161e-01;
            ++i;
            ref_result_100_hubner[i] = -1.865250e-01;
            ++i;
            ref_result_100_hubner[i] = -1.832132e-01;
            ++i;
            ref_result_100_hubner[i] = -1.797930e-01;
            ++i;
            ref_result_100_hubner[i] = -1.762763e-01;
            ++i;
            ref_result_100_hubner[i] = -1.726742e-01;
            ++i;
            ref_result_100_hubner[i] = -1.689971e-01;
            ++i;
            ref_result_100_hubner[i] = -1.652550e-01;
            ++i;
            ref_result_100_hubner[i] = -1.614569e-01;
            ++i;
            ref_result_100_hubner[i] = -1.576113e-01;
            ++i;
            ref_result_100_hubner[i] = -1.537257e-01;
            ++i;
            ref_result_100_hubner[i] = -1.498071e-01;
            ++i;
            ref_result_100_hubner[i] = -1.458615e-01;
            ++i;
            ref_result_100_hubner[i] = -1.418945e-01;
            ++i;
            ref_result_100_hubner[i] = -1.379106e-01;
            ++i;
            ref_result_100_hubner[i] = -1.339137e-01;
            ++i;
            ref_result_100_hubner[i] = -1.299069e-01;
            ++i;
            ref_result_100_hubner[i] = -1.258928e-01;
            ++i;
            ref_result_100_hubner[i] = -1.218729e-01;
            ++i;
            ref_result_100_hubner[i] = -1.178485e-01;
            ++i;
            ref_result_100_hubner[i] = -1.138196e-01;
            ++i;
            ref_result_100_hubner[i] = -1.097862e-01;
            ++i;
            ref_result_100_hubner[i] = -1.057471e-01;
            ++i;
            ref_result_100_hubner[i] = -1.017007e-01;
            ++i;
            ref_result_100_hubner[i] = -9.764491e-02;
            ++i;
            ref_result_100_hubner[i] = -9.357677e-02;
            ++i;
            ref_result_100_hubner[i] = -8.949286e-02;
            ++i;
            ref_result_100_hubner[i] = -8.538917e-02;
            ++i;
            ref_result_100_hubner[i] = -8.126107e-02;
            ++i;
            ref_result_100_hubner[i] = -7.710338e-02;
            ++i;
            ref_result_100_hubner[i] = -7.291035e-02;
            ++i;
            ref_result_100_hubner[i] = -6.867564e-02;
            ++i;
            ref_result_100_hubner[i] = -6.439234e-02;
            ++i;
            ref_result_100_hubner[i] = -6.005299e-02;
            ++i;
            ref_result_100_hubner[i] = -5.564948e-02;
            ++i;
            ref_result_100_hubner[i] = -5.117313e-02;
            ++i;
            ref_result_100_hubner[i] = -4.661461e-02;
            ++i;
            ref_result_100_hubner[i] = -4.196398e-02;
            ++i;
            ref_result_100_hubner[i] = -3.721056e-02;
            ++i;
            ref_result_100_hubner[i] = -3.234301e-02;
            ++i;
            ref_result_100_hubner[i] = -2.734922e-02;
            ++i;
            ref_result_100_hubner[i] = -2.221629e-02;
            ++i;
            ref_result_100_hubner[i] = -1.693047e-02;
            ++i;
            ref_result_100_hubner[i] = -1.147712e-02;
            ++i;
            ref_result_100_hubner[i] = -5.840619e-03;
            ++i;
            ref_result_100_hubner[i] = 0.000000e+00;

            DenseVector<DataType_> diff_2(test_line.copy());

            for(unsigned long i(0); i < g_w; ++i)
            {
                diff_2[i]= ref_result_100_hubner[i];
            }

            Difference<>::value(diff, test_line);
            Difference<>::value(diff_2, test_line);

            std::cout <<"Difference vector (GHIA): " << diff << std::endl;
            std::cout <<"Difference vector (HUEBNER): " << diff_2 << std::endl;

            double norm = Norm<vnt_l_two, false, Tag_>::value(diff);
            double norm_2 = Norm<vnt_l_two, false, Tag_>::value(diff_2);
            std::cout << "L2 norm (GHIA): " << norm << std::endl;
            std::cout << "L2 norm (HUEBNER): " << norm_2 << std::endl;
#ifdef DRIVEN_CAVITY_OUTPUT_ACCURACY
            std::string output_file = "accuracy_lbm_navsto_grid.dat";

            std::ofstream out_file_stream;

            out_file_stream.open(output_file.c_str(), ios::app);

            out_file_stream << "----------------------------------------------" << std::endl;
            out_file_stream << "timesteps = " << timesteps << std::endl;
            out_file_stream << "precision = " << sizeof(DataType_) << std::endl;
            out_file_stream << "  delta_x = " << dx << std::endl;
            out_file_stream << "  delta_y = " << dy << std::endl;
            out_file_stream << "  delta_t = " << dt << std::endl;
            out_file_stream << "      tau = " << tau << std::endl;
            out_file_stream << "       Re = " << _reynolds << std::endl;

            out_file_stream << "        U = " << lid_U << std::endl;

            out_file_stream << "L2 NORM(DIFF(GHIA)) = " << norm << std::endl;
            out_file_stream << "L2 NORM(DIFF(HUEBNER)) = " << norm_2 << std::endl;

            out_file_stream << "----------------------------------------------" << std::endl;;
            out_file_stream.close();
#endif

        }
};
//SolverLABNAVSTOGridDCTest<tags::CPU, double> solver_test_double("double", double(1), double(1), double(1.2), double(100), double(1), double(0));
//SolverLABNAVSTOGridDCTest<tags::CPU, float> solver_test_float("float", float(1), float(1), float(1.2), float(100), float(1), float(0));
//SolverLABNAVSTOGridDCTest<tags::CPU, double> solver_test_double_10("double", double(1.), double(1.), double(1.), double(100), double(1.), double(0));
