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
#include <honei/lbm/solver_labswe.hh>
#include <honei/swe/post_processing.hh>
#include <honei/swe/volume.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;

template <typename Tag_, typename DataType_>
class SolverLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_quick_test<" + type + ">")
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

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h, g_w, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> s_y(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver(1.,1.,1., g_w, g_h, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);
            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            std::cout << h << std::endl;
#endif
            TEST_CHECK(true);
        }

};
SolverLABSWETest<tags::CPU, float> solver_test_float("float");
SolverLABSWETest<tags::CPU, double> solver_test_double("double");
SolverLABSWETest<tags::CPU::MultiCore, float> solver_test_float_mc("float");
SolverLABSWETest<tags::CPU::MultiCore, double> solver_test_double_mc("double");
#ifdef HONEI_SSE
SolverLABSWETest<tags::CPU::SSE, float> solver_test_float_sse("float");
SolverLABSWETest<tags::CPU::SSE, double> solver_test_double_sse("double");
SolverLABSWETest<tags::CPU::MultiCore::SSE, float> solver_test_float_mc_sse("float");
SolverLABSWETest<tags::CPU::MultiCore::SSE, double> solver_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
SolverLABSWETest<tags::Cell, float> solver_test_float_cell("float");
SolverLABSWETest<tags::Cell, double> solver_test_double_cell("double");
#endif
#ifdef HONEI_CUDA
SolverLABSWETest<tags::GPU::CUDA, float> solver_test_float_cuda("float");
#endif
template <typename Tag_, typename DataType_>
class SolverLABSWEDrivenCavityTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEDrivenCavityTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_driven_cavity_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(129);
            unsigned long g_w(129);
            DataType_ p_d_t(1.);
            DataType_ p_d_x(1.);
            DataType_ p_d_y(1.);
            DataType_ veloc(0.129198966408);
            DataType_ tau(1.);
            unsigned long timesteps(1000);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            //set initial Dirichlet Boundaries:
            for(unsigned long i(0) ; i < g_w ; ++i)
            {
                u(0, i) = DataType_(veloc);
            }

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h, g_w, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> s_y(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::DRIVEN_CAVITY> solver(p_d_x, p_d_y, p_d_t, g_w, g_h, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);

            solver.set_relaxation_time(tau);
            solver.set_lid_velocity(veloc);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(u, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            std::cout << u << std::endl;
#endif
            TEST_CHECK(true);

            double reynolds(veloc * g_w /(p_d_x * p_d_x * (2. * tau - 1.)/ (6. * p_d_t)));
            std::cout<<"Re: " << reynolds << std::endl;

            DenseVector<DataType_> test_line(g_h);
            std::cout<<"Index: " << (g_w-1)/2 << std::endl;
            for(unsigned long i(0); i < g_h; ++i)
            {
                test_line[i] = DataType_(u(i,(g_w-1)/2)/veloc);
            }

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
                TEST_CHECK_EQUAL_WITHIN_EPS(test_line[indices_100[i]], ref_result_100[i], 0.5);
            }

            Difference<>::value(diff, test_line);

            std::cout <<"Difference vector: " << diff << std::endl;

        }

};

SolverLABSWEDrivenCavityTest<tags::CPU, double> solver_DV_test_double("double");

template <typename Tag_, typename DataType_>
class SolverLABSWEMassConservationTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEMassConservationTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(32);
            unsigned long g_w(32);
            unsigned long timesteps(5000);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cuboid<DataType_> a(h, DataType_(5.), DataType_(5.), DataType_(0.02), 25, 25);
            a.value();

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h, g_w, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> s_y(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver(1.,1.,1., g_w, g_h, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);
            solver.do_preprocessing();

#ifdef SOLVER_POSTPROCESSING_VOLUME
            DataType_ volumes[timesteps + 1];
            volumes[0] = 0.;
#endif

            for(unsigned long i(1); i <= timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif

#ifdef SOLVER_POSTPROCESSING_VOLUME
                volumes[i] = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(h, DataType_(0), DataType_(g_w), DataType_(1.), DataType_(1.));
#endif

            }
#ifdef SOLVER_VERBOSE
            std::cout << h << std::endl;
#endif

            DataType_ ana_vol(5. * 5. * 0.02 + g_w * g_h * 0.05);
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(h, DataType_(0), DataType_(g_w), DataType_(1.), DataType_(1.));
            std::cout << "Vol.: " << vol << std::endl;

#ifdef SOLVER_POSTPROCESSING_VOLUME
            std::string filename;
            std::ofstream file;
            filename = "out_labswe_vol.dat";
            file.open(filename.c_str());
            for(unsigned long i(1); i <= timesteps; ++i)
            {
                file << stringify(i) + " " + stringify(volumes[i]) + " " + stringify(fabs(volumes[i] - ana_vol))  +"\n";
            }
            file.close();
#endif

            TEST_CHECK_EQUAL_WITHIN_EPS(vol, ana_vol, 0.1);
        }

};
SolverLABSWEMassConservationTest<tags::CPU, double> solver_mass_cons_test_double("mass conservation, double");
SolverLABSWEMassConservationTest<tags::CPU, float> solver_mass_cons_test_float("mass conservation, float");
#ifdef HONEI_SSE
SolverLABSWEMassConservationTest<tags::CPU::SSE, double> solver_mass_cons_test_double_sse("mass conservation, double");
SolverLABSWEMassConservationTest<tags::CPU::SSE, float> solver_mass_cons_test_float_sse("mass conservation, float");
#endif
#ifdef HONEI_CELL
SolverLABSWEMassConservationTest<tags::Cell, double> solver_mass_cons_test_double_cell("mass conservation, double");
SolverLABSWEMassConservationTest<tags::Cell, float> solver_mass_cons_test_float_cell("mass conservation, float");
#endif
#ifdef HONEI_CUDA
SolverLABSWEMassConservationTest<tags::GPU::CUDA, float> solver_mass_cons_test_float_cuda("mass conservation, float");
#endif
template <typename Tag_, typename DataType_>
class SolverLABSWENOSLIPTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWENOSLIPTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswenoslip_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(1);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h, g_w, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> s_y(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(1.,1.,1., g_w, g_h, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);
            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            std::cout << h << std::endl;
#endif
            TEST_CHECK(true);
        }

};
SolverLABSWENOSLIPTest<tags::CPU, float> solver_noslip_test_float("float");
SolverLABSWENOSLIPTest<tags::CPU, double> solver_noslip_test_double("double");

