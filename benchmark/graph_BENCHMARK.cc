/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
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

#include <honei/graph/position.hh>
#include <honei/graph/graph.hh>
#include <honei/util/tags.hh>
#include <honei/graph/test_scenario.hh>

#include <string>

using namespace honei;

#define RUNS 3;
#define ITERATIONS 400;
template <typename Tag_, typename DataType_, typename GraphTag_>
class GraphBench :
    public Benchmark
{
    private:
    public:
        GraphBench(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            for (int size(1); size <= 30; ++size)
            {
            for(unsigned long i(0); i < 5; ++i)
            {
                Graph<DataType_> * g(TestScenario<DataType_>::Grid(size, size));
                Positions<Tag_, DataType_, GraphTag_> position(*g, DataType_(1));
                BENCHMARK(position.update(0.0, 400));
                
             //   std::cout <<"Run " << i<< "\tMaximum node force "<< position.max_node_force() << " after " << position.number_of_iterations() << " iterations.\n";
                delete g;
           }
                calculate();
                std::cout <<  size*size <<" " << _median << "\n";
                _benchlist.clear();
           }
        }
};



//GraphBench<tags::CPU, float, methods::WeightedFruchtermanReingold> wfr_graph_bench_float1("WFR Grid float BENCHMARK");
//GraphBench<tags::CPU::MultiCore, float, methods::WeightedFruchtermanReingold> wfr_graph_bench_float2("MC Grid float BENCHMARK");

#ifdef HONEI_SSE
//GraphBench<tags::CPU::SSE, float, methods::WeightedFruchtermanReingold> wfr_graph_bench_float3("SSE WFR Grid float BENCHMARK");
//GraphBench<tags::CPU::MultiCore::SSE, float, methods::WeightedFruchtermanReingold> wfr_graph_bench_float4("MC SSE WFR Grid float BENCHMARK");
#endif
#ifdef HONEI_CELL
//GraphBench<tags::Cell, float, methods::WeightedFruchtermanReingold> wfr_graph_bench_float5("Cell WFR Grid float BENCHMARK");
#endif
