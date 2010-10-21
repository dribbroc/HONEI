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

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

#include <string>
#endif

#include <honei/swe/relax_solver.hh>


using namespace std;
using namespace honei;

template <typename Tag_>
class RelaxSolverMIXEDPRECINNERBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        RelaxSolverMIXEDPRECINNERBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            ulint dwidth = _size;
            ulint dheight = _size;
            DenseMatrix<float> height(dheight, dwidth, float(5));
            for(ulint i = 0; i< height.rows(); ++i)
            {
                for(ulint j=height.columns()-10; j<height.columns(); ++j)
                {
                    height[i][j] = float(10);
                }
            }

            DenseMatrix<float> bottom(dheight, dwidth, float(1));

            DenseMatrix<float> u1(dheight, dwidth, float(0));
            DenseMatrix<float> u2(dheight, dwidth, float(0));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<float> u(entries, float(1));
            DenseVector<float> v(entries, float(1));
            DenseVector<float> w(entries, float(1)); 
            DenseVector<float> bx (entries/3, float(0));
            DenseVector<float> by (entries/3, float(0));
            DenseVector<float> c (3,float(5));
            DenseVector<float> d (3,float(5));
            c[0] = 10;
            c[1] = 6;
            c[2] = 11;
            d[0] = 10;
            d[1] = 5;
            d[2] = 11;

            float deltax = 5;
            float deltay = 5;
            float deltat = 5./24.;

            double eps = 10e-6;
            float manning = float(0);
            Scenario<float, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);
            scenario.height = &height;
            scenario.bottom = &bottom;
            scenario.x_veloc = &u1;
            scenario.y_veloc = &u2;
            scenario.u = &u;
            scenario.v = &v;
            scenario.w = &w;
            scenario.bottom_slopes_x = &bx;
            scenario.bottom_slopes_y = &by;
            scenario.c = &c;
            scenario.d = &d;
            scenario.delta_x = deltax;
            scenario.delta_y = deltay;
            scenario.delta_t = deltat;
            scenario.eps = eps;
            scenario.manning_n = manning;

            DenseVector<double> bx_2(bx.size());
            DenseVector<double> by_2(by.size());
            convert(bx_2, bx);
            convert(by_2, by);
            RelaxSolver<Tag_, float, float, double, float, float, source_types::SIMPLE, boundaries::REFLECT, MIXED> relax_solver
                (scenario, bx_2, by_2);
            relax_solver.do_preprocessing();
            for (ulint i = 1; i <= _count; ++i)
            {
                BENCHMARK(relax_solver.solve());
            }
            evaluate();

        }
};

RelaxSolverMIXEDPRECINNERBench<tags::CPU::MultiCore> mc_solver_bench_mp2("MC RelaxSolver Benchmark - size: 41, mixedprec 2", 41, 100);
RelaxSolverMIXEDPRECINNERBench<tags::CPU> solver_bench_mp22("RelaxSolver Benchmark - size: 41, mixedprec 2", 41, 100);
#ifdef HONEI_SSE
RelaxSolverMIXEDPRECINNERBench<tags::CPU::MultiCore::SSE> mc_sse_solver_bench_mp2("MC SSE RelaxSolver Benchmark - size: 41, mixedprec 2", 41, 100);
RelaxSolverMIXEDPRECINNERBench<tags::CPU::SSE> sse_solver_bench_mp2("SSE RelaxSolver Benchmark - size: 41, mixedprec 2", 41, 100);
#endif

