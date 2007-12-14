
/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
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

#include <libgraph/position.hh>
#include <unittest/unittest.hh>
#include <libla/dense_matrix.hh>
#include <libutil/tags.hh>

#include <string>


using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class KamadaKawaiPositionsBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        KamadaKawaiPositionsBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
        }

        virtual void run()
        {
            // Creatoing test scenario
            DataType_ pos[2*_nodecount];
            for (int i(0); i < _nodecount; ++i)
            {
                    pos[2*i+1] = sin((DataType_)i /(DataType_)_nodecount * 2.0f * 3.14f);
                    pos[2*i] = cos((DataType_)i / (DataType_)_nodecount * 2.0f * 3.14f);
            }

            DataType_ adj[_nodecount*_nodecount];
            for (int i(0); i < _nodecount; ++i)
                for (int j(0); j < _nodecount; ++j)
                    adj[i*_nodecount + j] = i == j ? 0 : 1;

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(_nodecount,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<bool> > pNeighbour(new SparseMatrix<bool>(_nodecount,_nodecount));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::KamadaKawai> position(*pPosition, *pNeighbour, 2);

            // update the positions
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(position.update(0.00001, _nodecount * 10));
            }

            std::cout << "max_node_force of KK  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of KK  "<< position.number_of_iterations() << std::endl;

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class FruchtermanReingoldPositionsBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        FruchtermanReingoldPositionsBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
        }

        virtual void run()
        {
            // Creatoing test scenario
            DataType_ pos[2*_nodecount];
            for (int i(0); i < _nodecount; ++i)
            {
                    pos[2*i+1] = sin((DataType_)i /(DataType_)_nodecount * 2.0f * 3.14f);
                    pos[2*i] = cos((DataType_)i / (DataType_)_nodecount * 2.0f * 3.14f);
            }

            DataType_ adj[_nodecount*_nodecount];
            for (int i(0); i < _nodecount; ++i)
                for (int j(0); j < _nodecount; ++j)
                    adj[i*_nodecount + j] = i == j ? 0 : 1;

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(_nodecount,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<bool> > pNeighbour(new SparseMatrix<bool>(_nodecount,_nodecount));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::FruchtermanReingold> position(*pPosition, *pNeighbour, 2);

            // update the positions
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(position.update(0.00001, _nodecount * 5));
            }

            std::cout << "max_node_force of FR  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of FR  "<< position.number_of_iterations() << std::endl;

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class WeightedKamadaKawaiPositionsBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        WeightedKamadaKawaiPositionsBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
        }

        virtual void run()
        {
            // Creatoing test scenario
            DataType_ pos[2*_nodecount];
            for (int i(0); i < _nodecount; ++i)
            {
                    pos[2*i+1] = sin((DataType_)i /(DataType_)_nodecount * 2.0f * 3.14f);
                    pos[2*i] = cos((DataType_)i / (DataType_)_nodecount * 2.0f * 3.14f);
            }

            DataType_ edge_weights[_nodecount*_nodecount];
            for (int i(0); i < _nodecount; ++i)
                for (int j(0); j < _nodecount; ++j)
                    edge_weights[i*_nodecount + j] = i == j ? 0 : 1;

            DataType_ node_weights[_nodecount];
            for (int i(0); i < _nodecount; ++i)
            {
                    node_weights[i] = 1;
            }

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(_nodecount,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++]; 
            }

            i = 0;
            std::tr1::shared_ptr<DenseVector<DataType_> > pNode_Weights(new DenseVector<DataType_>(_nodecount));
            for (typename Vector<DataType_>::ElementIterator e(pNode_Weights->begin_elements()),
                    e_end(pNode_Weights->end_elements()); e != e_end ; ++e)
            {
                *e = 3*node_weights[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pEdge_Weights(new SparseMatrix<DataType_>(_nodecount,_nodecount));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEdge_Weights->begin_elements()),
                    e_end(pEdge_Weights->end_elements()); e != e_end ; ++e)
            {
                if (edge_weights[i] > std::numeric_limits<DataType_>::epsilon()) *e = edge_weights[i];
                i++;
            }

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, *pNode_Weights, *pEdge_Weights);

            // update the positions
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(position.update(0.00001, _nodecount * 10));
            }

            std::cout << "max_node_force of WKK  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of WKK  "<< position.number_of_iterations() << std::endl;

            evaluate();
        }
};

KamadaKawaiPositionsBench<tags::CPU, float> kamada_kawai_positions_bench_float("float", 200, 10);
KamadaKawaiPositionsBench<tags::CPU, double> kamada_kawai_positions_bench_double("double", 200, 10);
FruchtermanReingoldPositionsBench<tags::CPU, float> fruchterman_reingold_positions_bench_float("float", 200, 10);
FruchtermanReingoldPositionsBench<tags::CPU, double> fruchterman_reingold_positions_bench_double("double", 200, 10);
WeightedKamadaKawaiPositionsBench<tags::CPU, float> weighted_kamada_kawai_positions_bench_float("float", 200, 10);
WeightedKamadaKawaiPositionsBench<tags::CPU, double> weighted_kamada_kawai_positions_bench_double("double", 200, 10);
#ifdef HONEI_SSE
KamadaKawaiPositionsBench<tags::CPU::SSE, float> sse_kamada_kawai_positions_bench_float("SSE float", 200, 10);
KamadaKawaiPositionsBench<tags::CPU::SSE, double> sse_kamada_kawai_positions_bench_double("SSE double", 200, 10);
FruchtermanReingoldPositionsBench<tags::CPU::SSE, float> sse_fruchterman_reingold_positions_bench_float("SSE float", 200, 10);
FruchtermanReingoldPositionsBench<tags::CPU::SSE, double> sse_fruchterman_reingold_positions_bench_double("SSE double", 200, 10);
WeightedKamadaKawaiPositionsBench<tags::CPU::SSE, float> sse_weighted_kamada_kawai_positions_bench_float("SSE float", 200, 10);
WeightedKamadaKawaiPositionsBench<tags::CPU::SSE, double> sse_weighted_kamada_kawai_positions_bench_double("SSE double", 200, 10);
#endif
