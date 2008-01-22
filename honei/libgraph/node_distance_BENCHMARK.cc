/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 *
 * This file is part of the Graph C library. LibGraph is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
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

#include <honei/libgraph/node_distance.hh>


using namespace std;
using namespace honei;

namespace NodeDistanceMethods
{
    struct ForKamadaKawai
    {
    };

    struct ForFruchtermanReingold
    {
    };

    struct ForWeightedFruchtermanReingold
    {
    };
}

template <typename Tag_, typename DataType_, typename NDM_>
class NodeDistanceBench;

template <typename Tag_, typename DataType_>
class NodeDistanceBench<Tag_, DataType_, NodeDistanceMethods::ForKamadaKawai> :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        NodeDistanceBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a test scenario
            DenseMatrix<DataType_>  pPos(_nodecount, 2);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPos.begin_elements()), e_end(pPos.end_elements()); e != e_end ; ++e)
            {
                *e = e.index();
            }

            // Calculate distance matrix
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(NodeDistance<Tag_>::value(pPos));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class NodeDistanceBench<Tag_, DataType_, NodeDistanceMethods::ForFruchtermanReingold> :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        NodeDistanceBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a test scenario
            DenseMatrix<DataType_>  pPos(_nodecount, 2);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPos.begin_elements()), e_end(pPos.end_elements()); e != e_end ; ++e)
            {
                *e = e.index();
            }

            DenseMatrix<bool>  pAdj(_nodecount, _nodecount, false);
            for (typename MutableMatrix<bool>::ElementIterator e(pAdj.begin_elements()), e_end(pAdj.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column())
                {
                    *e = true;
                }
            }

            DenseMatrix<DataType_> sd(_nodecount, _nodecount);
            DenseMatrix<DataType_> isd(_nodecount, _nodecount);
            DataType_ rfr(_nodecount * 3);

            // Calculate distance matrix
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(NodeDistance<Tag_>::value(pPos, pAdj, sd, isd, rfr));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class NodeDistanceBench<Tag_, DataType_, NodeDistanceMethods::ForWeightedFruchtermanReingold> :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        NodeDistanceBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a test scenario
            DenseMatrix<DataType_>  pPos(_nodecount, 2);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPos.begin_elements()), e_end(pPos.end_elements()); e != e_end ; ++e)
            {
                *e = e.index();
            }

            DenseMatrix<DataType_>  pEW(_nodecount, _nodecount, DataType_(0));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEW.begin_elements()), e_end(pEW.end_elements()); e != e_end ; ++e)
            {
                if (e.row() < e.column())
                {
                    *e = DataType_(2 + (e.index() % 2));
                    pEW[e.column()][e.row()] = *e;
                }
            }

            DenseMatrix<DataType_> sd(_nodecount, _nodecount);
            DenseMatrix<DataType_> isd(_nodecount, _nodecount);
            DataType_ rfr(_nodecount * 3);

            // Calculate distance matrix
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(NodeDistance<Tag_>::value(pPos, pEW, sd, isd, rfr));
            }

            evaluate();
        }
};

NodeDistanceBench<tags::CPU, float, NodeDistanceMethods::ForKamadaKawai> node_distance_bench_float_KK("NodeDistance for KK Benchmark float", 1000, 10);
NodeDistanceBench<tags::CPU, double, NodeDistanceMethods::ForKamadaKawai> node_distance_bench_double_KK("NodeDistance for KK Benchmark double", 1000, 10);
NodeDistanceBench<tags::CPU, float, NodeDistanceMethods::ForFruchtermanReingold> node_distance_bench_float_FR("NodeDistance for FR Benchmark float", 1000, 10);
NodeDistanceBench<tags::CPU, double, NodeDistanceMethods::ForFruchtermanReingold> node_distance_bench_double_FR("NodeDistance for FR Benchmark double", 1000, 10);
NodeDistanceBench<tags::CPU, float, NodeDistanceMethods::ForWeightedFruchtermanReingold> node_distance_bench_float_WFR("NodeDistance for WFR Benchmark float", 0100, 10);
NodeDistanceBench<tags::CPU, double, NodeDistanceMethods::ForWeightedFruchtermanReingold> node_distance_bench_double_WFR("NodeDistance for WFR Benchmark double", 1000, 10);

#ifdef HONEI_SSE
NodeDistanceBench<tags::CPU::SSE, float, NodeDistanceMethods::ForKamadaKawai> sse_node_distance_bench_float_KK("SSE NodeDistance for KK Benchmark float", 1000, 10);
NodeDistanceBench<tags::CPU::SSE, double, NodeDistanceMethods::ForKamadaKawai> sse_node_distance_bench_double_KK("SSE NodeDistance for KK Benchmark double", 1000, 10);
NodeDistanceBench<tags::CPU::SSE, float, NodeDistanceMethods::ForFruchtermanReingold> sse_node_distance_bench_float_FR("SSE NodeDistance for FR Benchmark float", 1000, 10);
NodeDistanceBench<tags::CPU::SSE, double, NodeDistanceMethods::ForFruchtermanReingold> sse_node_distance_bench_double_FR("SSE NodeDistance for FR Benchmark double", 1000, 10);
NodeDistanceBench<tags::CPU::SSE, float, NodeDistanceMethods::ForWeightedFruchtermanReingold> sse_node_distance_bench_float_WFR("SSE NodeDistance for WFR Benchmark float", 1000, 10);
NodeDistanceBench<tags::CPU::SSE, double, NodeDistanceMethods::ForWeightedFruchtermanReingold> sse_node_distance_bench_double_WFR("SSE NodeDistance for WFR Benchmark double", 1000, 10);
#endif

#ifdef HONEI_CELL
NodeDistanceBench<tags::Cell, float, NodeDistanceMethods::ForKamadaKawai> cell_node_distance_bench_float_KK("Cell NodeDistance for KK Benchmark float", 1000, 10);
NodeDistanceBench<tags::Cell, float, NodeDistanceMethods::ForFruchtermanReingold> cell_node_distance_bench_float_FR("Cell NodeDistance for FR Benchmark float", 1000, 10);
NodeDistanceBench<tags::Cell, float, NodeDistanceMethods::ForWeightedFruchtermanReingold> cell_node_distance_bench_float_WFR("Cell NodeDistance for WFR Benchmark float", 1000, 10);
#endif