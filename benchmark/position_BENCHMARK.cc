
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

#include <honei/graph/position.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/util/tags.hh>

#include <string>


using namespace std;
using namespace honei;

namespace Scenarios
{
    struct Clique
    {
        static std::string get_name()
        {
            std::string name = "Clique";
            return name;
        }
    };

    struct SquareGrid
    {
        static std::string get_name()
        {
            std::string name = "SquareGrid";
            return name;
        }
    };

    struct BinaryTree
    {
        static std::string get_name()
        {
            std::string name = "BinaryTree";
            return name;
        }
    };
}

template <typename DataType_, typename ST_>
struct Scenario;

template <typename DataType_>
struct Scenario<DataType_, Scenarios::Clique>
{
    static void create(DenseMatrix<DataType_> & Position, DenseVector<DataType_> & Node_Weights, 
    SparseMatrix<DataType_> & Edge_Weights)
    {
        unsigned long _nodecount(Position.rows());

        // create Position
        for (typename DenseMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
                    e_end(Position.end_elements());e != e_end ; ++e)
            {
                e.column() == 0 ? *e = cos((DataType_) (e.row()) / (DataType_)_nodecount * 2.0f * 3.14f) :
                *e = sin((DataType_) (e.row()) /(DataType_)_nodecount * 2.0f * 3.14f);
            }

        // create Node_Weights
        for (typename Vector<DataType_>::ElementIterator e(Node_Weights.begin_elements()),
                e_end(Node_Weights.end_elements()); e != e_end ; ++e)
        {
            *e = 1;
        }

        // create Edge_Weights
        for (typename DenseMatrix<DataType_>::ElementIterator e(Edge_Weights.begin_elements()),
                    e_end(Edge_Weights.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column()) *e = 1;
            }
    }
};

template <typename DataType_>
struct Scenario<DataType_, Scenarios::SquareGrid>
{
    static void create(DenseMatrix<DataType_> & Position, DenseVector<DataType_> & Node_Weights, 
    SparseMatrix<DataType_> & Edge_Weights)
    {
        unsigned long _nodecount(Position.rows());
        unsigned long _nodecount_2((unsigned long)sqrt(_nodecount));

        // create Position
        for (typename DenseMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
                    e_end(Position.end_elements());e != e_end ; ++e)
            {
                e.column() == 0 ? *e = cos((DataType_) (e.row()) / (DataType_)_nodecount * 2.0f * 3.14f) :
                *e = sin((DataType_) (e.row()) /(DataType_)_nodecount * 2.0f * 3.14f);
            }

        // create Node_Weights
        for (typename Vector<DataType_>::ElementIterator e(Node_Weights.begin_elements()),
                e_end(Node_Weights.end_elements()); e != e_end ; ++e)
        {
            *e = 1;
        }

        // create Edge_Weights
        for (typename MutableMatrix<DataType_>::ElementIterator e(Edge_Weights.begin_elements()),
                    e_end(Edge_Weights.end_elements()); e != e_end ; ++e)
            {
                if ((((e.row() + 1 == e.column()) && (e.column() % _nodecount_2)) || (e.column() == e.row() + _nodecount_2)) && 
                (e.row() != e.column())) 
                {
                    *e = 1;
                    Edge_Weights[e.column()][e.row()] = 1;
                }
            }
    }
};

template <typename DataType_>
struct Scenario<DataType_, Scenarios::BinaryTree>
{
    static void create(DenseMatrix<DataType_> & Position, DenseVector<DataType_> & Node_Weights, 
    SparseMatrix<DataType_> & Edge_Weights)
    {
        unsigned long _nodecount(Position.rows());

        // create Position
        for (typename DenseMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
                    e_end(Position.end_elements());e != e_end ; ++e)
            {
                e.column() == 0 ? *e = cos((DataType_) (e.row()) / (DataType_)_nodecount * 2.0f * 3.14f) :
                *e = sin((DataType_) (e.row()) /(DataType_)_nodecount * 2.0f * 3.14f);
            }

        // create Node_Weights
        for (typename Vector<DataType_>::ElementIterator e(Node_Weights.begin_elements()),
                e_end(Node_Weights.end_elements()); e != e_end ; ++e)
        {
            *e = 1;
        }

        // create Edge_Weights
        for (typename DenseMatrix<DataType_>::ElementIterator e(Edge_Weights.begin_elements()),
                    e_end(Edge_Weights.end_elements()); e != e_end ; ++e)
            {
                if ((e.column() == e.row() * 2 +1) || (e.column() == e.row() * 2 +2)) 
                {
                    *e = 1;
                    Edge_Weights[e.column()][e.row()] = 1;
                }
            }
    }
};

template <typename Tag_, typename DataType_, typename ST_>
class WeightedKamadaKawaiPositionsBench :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        WeightedKamadaKawaiPositionsBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            if (ST_::get_name() == "Clique") _nodecount = nodecount;
            if (ST_::get_name() == "SquareGrid") _nodecount = nodecount*nodecount;
            if (ST_::get_name() == "BinaryTree") _nodecount = (1 << (nodecount + 1)) -1;
            _count = count;
        }

        virtual void run()
        {
            // Creatoing test scenario
            DenseMatrix<DataType_> Coordinates(_nodecount, 2, DataType_(0));
            DenseVector<DataType_> Node_Weights(_nodecount, DataType_(0));
            SparseMatrix<DataType_> Edge_Weights(_nodecount, _nodecount);
            Scenario<DataType_, ST_>::create(Coordinates, Node_Weights, Edge_Weights);

            // Creating a Positions object with the test scenario and update the positions
            unsigned long number_of_iterations;
            DataType_ max_node_force;
            DenseMatrix<DataType_> new_coordinates(_nodecount, 2, DataType_(0));
            DenseMatrix<DataType_> pos_copy(Coordinates.copy());
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(pos_copy, Node_Weights, Edge_Weights);
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(position.update(0, _nodecount * 10));
                number_of_iterations = position.number_of_iterations();
                max_node_force = position.max_node_force();
                new_coordinates = position.coordinates();
                pos_copy = Coordinates.copy();
            }
            std::cout << "number of nodes of WKK  "<< _nodecount << std::endl;
            std::cout << "Maximum node force of WKK:  "<< max_node_force << std::endl;
            std::cout << "Number of iterations of WKK:  "<< number_of_iterations << std::endl;

            evaluate();
        }
};

template <typename Tag_, typename DataType_, typename ST_>
class WeightedFruchtermanReingoldPositionsBench :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        WeightedFruchtermanReingoldPositionsBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            if (ST_::get_name() == "Clique") _nodecount = nodecount;
            if (ST_::get_name() == "SquareGrid") _nodecount = nodecount*nodecount;
            if (ST_::get_name() == "BinaryTree") _nodecount = (1 << (nodecount + 1)) -1;
            _count = count;
        }

        virtual void run()
        {
            // Creatoing test scenario
            DenseMatrix<DataType_> Coordinates(_nodecount, 2, DataType_(0));
            DenseVector<DataType_> Node_Weights(_nodecount, DataType_(0));
            SparseMatrix<DataType_> Edge_Weights(_nodecount, _nodecount);
            Scenario<DataType_, ST_>::create(Coordinates, Node_Weights, Edge_Weights);

            // Creating a Positions object with the test scenario and update the positions
            unsigned long number_of_iterations;
            DataType_ max_node_force;
            DenseMatrix<DataType_> new_coordinates(_nodecount, 2, DataType_(0));
            DenseMatrix<DataType_> pos_copy(Coordinates.copy());
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(pos_copy, Node_Weights, Edge_Weights);                
            for(unsigned long i = 0; i < _count; ++i)
            {                
                BENCHMARK(position.update(0, _nodecount * 5));
                number_of_iterations = position.number_of_iterations();
                max_node_force = position.max_node_force();
                new_coordinates = position.coordinates();
                pos_copy = Coordinates.copy();
            }
            std::cout << "number of nodes of WFR  "<< _nodecount << std::endl;
            std::cout << "Maximum node force of WFR:  "<< max_node_force << std::endl;
            std::cout << "Number of iterations of WFR:  "<< number_of_iterations << std::endl;

            evaluate();
        }
};

#define POSITIONBENCH Scenarios::SquareGrid // possible scenarios are: Clique, SquareGrid, BinaryTree
#define POSITIONBENCHSIZE 20 //POSITIONBENCHSIZE = numbers of nodes (Clique), POSITIONBENCHSIZE = numbers of nodes in a line (SquareGrid), POSITIONBENCHSIZE = depth (BinaryTree)
#define POSITIONBENCHCOUNT 3

WeightedKamadaKawaiPositionsBench<tags::CPU, float, POSITIONBENCH> weighted_kamada_kawai_positions_bench_float("WeightedKamadaKawai Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU, double, POSITIONBENCH> weighted_kamada_kawai_positions_bench_double("WeightedKamadaKawai Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU, float, POSITIONBENCH> weighted_fruchterman_reingold_positions_bench_float("WeightedFruchtermanReingold Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU, double, POSITIONBENCH> weighted_fruchterman_reingold_positions_bench_double("WeightedFruchtermanReingold Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU::MultiCore, float, POSITIONBENCH> mc_weighted_kamada_kawai_positions_bench_float("MC WeightedKamadaKawai Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU::MultiCore, double, POSITIONBENCH> mc_weighted_kamada_kawai_positions_bench_double("MC WeightedKamadaKawai Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::MultiCore, float, POSITIONBENCH> mc_weighted_fruchterman_reingold_positions_bench_float("MC WeightedFruchtermanReingold Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::MultiCore, double, POSITIONBENCH> mc_weighted_fruchterman_reingold_positions_bench_double("MC WeightedFruchtermanReingold Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);

#ifdef HONEI_SSE
WeightedKamadaKawaiPositionsBench<tags::CPU::SSE, float, POSITIONBENCH> sse_weighted_kamada_kawai_positions_bench_float("SSE WeightedKamadaKawai Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU::SSE, double, POSITIONBENCH> sse_weighted_kamada_kawai_positions_bench_double("SSE WeightedKamadaKawai Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::SSE, float, POSITIONBENCH> sse_weighted_fruchterman_reingold_positions_bench_float("SSE WeightedFruchtermanReingold Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::SSE, double, POSITIONBENCH> sse_weighted_fruchterman_reingold_positions_bench_double("SSE WeightedFruchtermanReingold Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU::MultiCore::SSE, float, POSITIONBENCH> mc_sse_weighted_kamada_kawai_positions_bench_float("MC SSE WeightedKamadaKawai Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedKamadaKawaiPositionsBench<tags::CPU::MultiCore::SSE, double, POSITIONBENCH> mc_sse_weighted_kamada_kawai_positions_bench_double("MC SSE WeightedKamadaKawai Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::MultiCore::SSE, float, POSITIONBENCH> mc_sse_weighted_fruchterman_reingold_positions_bench_float("MC SSE WeightedFruchtermanReingold Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::CPU::MultiCore::SSE, double, POSITIONBENCH> mc_sse_weighted_fruchterman_reingold_positions_bench_double("MC SSE WeightedFruchtermanReingold Benchmark double", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
#endif
#ifdef HONEI_CELL
WeightedKamadaKawaiPositionsBench<tags::Cell, float, POSITIONBENCH> cell_weighted_kamada_kawai_positions_bench_float("Cell WeightedKamadaKawai Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
WeightedFruchtermanReingoldPositionsBench<tags::Cell, float, POSITIONBENCH> cell_weighted_fruchterman_reingold_positions_bench_float("Cell WeightedFruchtermanReingold Benchmark float", POSITIONBENCHSIZE, POSITIONBENCHCOUNT);
#endif
