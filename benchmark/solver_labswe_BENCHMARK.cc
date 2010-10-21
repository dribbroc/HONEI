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

#include <string>
#endif

#include <honei/lbm/solver_labswe.hh>
#include <honei/swe/volume.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class LBMSolverBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        LBMSolverBench(const std::string & id, unsigned long size, int count) :
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

            SolverLABSWE<Tag_, DataType_, lbm_force::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver(1.,1.,1., g_w, g_h, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);
            solver.do_preprocessing();

            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(solver.solve());
            }
            evaluate();
        }
};

LBMSolverBench<tags::CPU, float> solver_bench_float_1("LBM solver Benchmark - size: 50, float", 50, 100);
LBMSolverBench<tags::CPU, double> solver_bench_double_1("LBM solver Benchmark - size: 50, double", 50, 100);
#ifdef HONEI_SSE
LBMSolverBench<tags::CPU::SSE, float> solver_bench_float_sse_1("LBM solver Benchmark - size: 50, float, SSE", 50, 100);
LBMSolverBench<tags::CPU::SSE, double> solver_bench_double_sse_1("LBM solver Benchmark - size: 50, double, SSE", 50, 100);
#endif
