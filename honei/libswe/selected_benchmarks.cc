/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <benchmark/benchmark.cc>

#include <honei/libswe/relax_solver.hh>


using namespace honei;


template <typename Tag_, typename DataType_>
class RelaxSolverBench :
    public Benchmark
{
    private:
        int _count;
    public:
        RelaxSolverBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _count = count;
        _plots = true;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long counter(10) ; counter <= _count ; counter+=10)
            {
                unsigned long _size(counter);
                cores.push_back(Tag_::name);
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
                /*for(ulint i = 0; i< bottom.rows(); ++i)
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
                  }*/

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
                DataType_ deltat = 5./24.;

                double eps = 10e-6;
                DataType_ manning = DataType_(0);
                Scenario<DataType_, swe_solvers::RELAX, boundaries::REFLECT> scenario(dwidth, dheight);
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

                RelaxSolver<Tag_, DataType_, DataType_, DataType_, DataType_, DataType_, source_types::SIMPLE, boundaries::REFLECT> relax_solver
                    (scenario);
                relax_solver.do_preprocessing();
                string outHeight = stringify(height);

                for(int i = 0; i < 30; ++i)
                {
                    BENCHMARK(relax_solver.solve());
                }
                info.size.clear();
                info.size.push_back(_size);
                infolist.push_back(info);
                std::cout << "finished run " << counter << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 30);
        }
};

#ifndef HONEI_CELL
RelaxSolverBench<tags::CPU, float> solver_bench_float_2("RelaxSolver Benchmark - float", 120);
RelaxSolverBench<tags::CPU, double> solver_bench_double_2("RelaxSolver Benchmark - double", 120);
RelaxSolverBench<tags::CPU::MultiCore, float> mc_solver_bench_float_2("MC RelaxSolver Benchmark - float", 120);
RelaxSolverBench<tags::CPU::MultiCore, double> mc_solver_bench_double_2("MC RelaxSolver Benchmark - double", 120);
#ifdef HONEI_SSE
RelaxSolverBench<tags::CPU::MultiCore::SSE, float> mc_sse_solver_bench_float_2("MC SSE RelaxSolver Benchmark - float", 120);
RelaxSolverBench<tags::CPU::MultiCore::SSE, double> mc_sse_solver_bench_double_2("MC SSE RelaxSolver Benchmark double", 120);
RelaxSolverBench<tags::CPU::SSE, float> sse_solver_bench_float_2("SSE RelaxSolver Benchmark - float", 120);
RelaxSolverBench<tags::CPU::SSE, double> sse_solver_bench_double_2("SSE RelaxSolver Benchmark - double", 120);
#endif
#elif HONEI_CELL
RelaxSolverBench<tags::Cell, float> cell_solver_bench_float_2("Cell RelaxSolverSum Benchmark - float",
         120);
#endif

