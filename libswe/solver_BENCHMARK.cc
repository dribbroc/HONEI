
/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <tr1/memory>
#include <string>
#endif

#include <libswe/solver.hh>


using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>

class RelaxSolverBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        RelaxSolverBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            ulint dwidth = _size;
            ulint dheight = _size;

            DenseMatrix<DataType_> height(dheight, dwidth, DataType_(5));
            //SCENARIO setup
            for(ulint i = 0; i< height.rows(); ++i)
            {
                for(ulint j=height.columns()-10; j<height.columns(); ++j)
                {
                    height[i][j] = DataType_(10);
                }
                 //(height)[0][i] = DataType_(10);
            }
            //END SCENARIO setup
            DenseMatrix<DataType_> bottom(dheight, dwidth, DataType_(1));
            for(ulint i = 0; i< bottom.rows(); ++i)
            {
                for(ulint j=0; j<bottom.columns()-10; ++j)
                {
                    bottom[i][j] = DataType_(1);
                    if(j>4 && j< bottom.columns()-9)
                    {
                        if(i < 6 || i > 11)
                            bottom[i][j] = DataType_(3);
                        else
                            bottom[i][j] = DataType_(1);
                    }
                }
            }

            DenseMatrix<DataType_> u1(dheight, dwidth, DataType_(0));
            DenseMatrix<DataType_> u2(dheight, dwidth, DataType_(0));
            unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
            DenseVector<DataType_> u(entries, DataType_(1));
            DenseVector<DataType_> v(entries, DataType_(1));
            DenseVector<DataType_> w(entries, DataType_(1)); 
            DenseVector<DataType_> bx (entries/3, DataType_(0));
            DenseVector<DataType_> by (entries/3, DataType_(0));
            DenseVector<DataType_> c (3,DataType_(5));
            DenseVector<DataType_> d (3,DataType_(5));
            //SCENARIO setup:
            c[0] = 10;
            c[1] = 6;
            c[2] = 11;
            d[0] = 10;
            d[1] = 5;
            d[2] = 11;

            DataType_ deltax = 5;
            DataType_ deltay = 5;
            DataType_ deltat = 5./22.;

            double eps = 10e-6;
            DataType_ manning = DataType_(0);

            RelaxSolver<Tag_, DataType_, DataType_, DataType_, DataType_, DataType_> relax_solver
                (&height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);
            relax_solver.do_preprocessing();

            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(relax_solver.solve());
            }
            //BenchmarkInfo info(Sum<>::get_benchmark_info(dv0, dv1));
            //evaluate(info);
            evaluate();
        }
};

RelaxSolverBench<tags::CPU, float> relax_solver_bench_float_1("RelaxSolver Benchmark - size: 21, float", 21, 100);
RelaxSolverBench<tags::CPU, double> relax_solver_bench_double_1("RelaxSolver Benchmark - size: 21, double", 21, 100);
RelaxSolverBench<tags::CPU, float> relax_solver_bench_float_2("RelaxSolver Benchmark - size: 41, float", 41, 100);
RelaxSolverBench<tags::CPU, double> relax_solver_bench_double_2("RelaxSolver Benchmark - size: 41, double", 41, 100);
RelaxSolverBench<tags::CPU::MultiCore, float> mc_relax_solver_bench_float_1("MC RelaxSolver Benchmark - size: 21, float", 21, 100);
RelaxSolverBench<tags::CPU::MultiCore, double> mc_relax_solver_bench_double_1("MC RelaxSolver Benchmark - size: 21, double", 21, 100);
RelaxSolverBench<tags::CPU::MultiCore, float> mc_relax_solver_bench_float_2("MC RelaxSolver Benchmark - size: 41, float", 41, 100);
RelaxSolverBench<tags::CPU::MultiCore, double> mc_relax_solver_bench_double_2("MC RelaxSolver Benchmark - size: 41, double", 41, 100);
#ifdef HONEI_SSE
RelaxSolverBench<tags::CPU::MultiCore::SSE, float> mc_sse_relax_solver_bench_float_1("MC SSE RelaxSolver Benchmark - size: 21, float", 21, 100);
RelaxSolverBench<tags::CPU::MultiCore::SSE, double> mc_sse_relax_solver_bench_double_1("MC SSE RelaxSolver Benchmark - size: 21, double", 21, 100);
RelaxSolverBench<tags::CPU::MultiCore::SSE, float> mc_sse_relax_solver_bench_float_2("MC SSE RelaxSolver Benchmark - size: 41, float", 41, 100);
RelaxSolverBench<tags::CPU::MultiCore::SSE, double> mc_sse_relax_solver_bench_double_2("MC SSE RelaxSolver Benchmark - size: 41, double", 41, 100);
RelaxSolverBench<tags::CPU::SSE, float> sse_relax_solver_bench_float_1("SSE RelaxSolver Benchmark - size: 21, float", 21, 100);
RelaxSolverBench<tags::CPU::SSE, double> sse_relax_solver_bench_double_1("SSE RelaxSolver Benchmark - size: 21, double", 21, 100);
RelaxSolverBench<tags::CPU::SSE, float> sse_relax_solver_bench_float_2("SSE RelaxSolver Benchmark - size: 41, float", 41, 100);
RelaxSolverBench<tags::CPU::SSE, double> sse_relax_solver_bench_double_2("SSE RelaxSolver Benchmark - size: 41, double", 41, 100);
#endif
#ifdef HONEI_CELL
RelaxSolverBench<tags::Cell, float> cell_relax_solver_bench_float_1("Cell RelaxSolverSum Benchmark - size: 21, float",
        21, 100);
RelaxSolverBench<tags::Cell, float> cell_relax_solver_bench_float_2("Cell RelaxSolverSum Benchmark - size: 41, float",
        41, 100);
#endif

