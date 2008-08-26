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

#include <honei/graph/breadth_first_search.hh>


using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class BreadthFirstSearchWeightedCliqueBench :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _count;
    public:
        BreadthFirstSearchWeightedCliqueBench(const std::string & id, unsigned long nodecount, unsigned long count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a nodecount-clique
            SparseMatrix<DataType_>  pEW(_nodecount, _nodecount);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEW.begin_elements()), e_end(pEW.end_elements()); e != e_end ; ++e)
            {
                if (e.row() < e.column())
                {
                    *e = DataType_(2 + (e.index() % 2));
                    pEW[e.column()][e.row()] = *e;
                }
            }

            DenseVector<DataType_>  pNW(_nodecount);
            for (typename DenseVector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(4);
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Calculate distance matrix
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class BreadthFirstSearchWeightedBinaryTreeBench :
    public Benchmark
{
    private:
        unsigned long _nodecount;
        unsigned long _depth;
        unsigned long _count;
    public:
        BreadthFirstSearchWeightedBinaryTreeBench(const std::string & id, unsigned long depth, unsigned long count) :
            Benchmark(id)
        {
            _count = count;
            _depth = depth;
            _nodecount = (1 << (_depth + 1)) -1;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a binary tree with deep _nodecount
            SparseMatrix<DataType_>  pEW(_nodecount, _nodecount);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                pEW[i][column] = DataType_(2);
                pEW[column][i] = DataType_(2);
                pEW[i][column + 1] = DataType_(2);
                pEW[column + 1][i] = DataType_(2);
            }

            DenseVector<DataType_>  pNW(_nodecount);
            unsigned long exp(0);
            for (typename DenseVector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(1 + (exp % 2));
                if (e.index() + 1 == (1 << (exp +1)) -1) exp++;
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Calculate distance matrix
            for(unsigned long i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW));
            }

            evaluate();
        }
};

BreadthFirstSearchWeightedCliqueBench<tags::CPU, float> breadth_first_search_weighted_clique_bench_float("BFS WeightedClique Benchmark float", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU, double> breadth_first_search_weighted_clique_bench_double("BFS WeightedClique Benchmark double", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU::MultiCore, float> mc_breadth_first_search_weighted_clique_bench_float("MC BFS WeightedClique Benchmark float", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU::MultiCore, double> mc_breadth_first_search_weighted_clique_bench_double("MC BFS WeightedClique Benchmark double", 1000, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU, float> breadth_first_search_weighted_binary_tree_bench_float("BFS WeightedBinaryTree Benchmark float", 10, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU, double> breadth_first_search_weighted_binary_tree_bench_double("BFS WeightedBinaryTree Benchmark double", 10, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU::MultiCore, float> mc_breadth_first_search_weighted_binary_tree_bench_float("MC BFS WeightedBinaryTree Benchmark float", 10, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU::MultiCore, double> mc_breadth_first_search_weighted_binary_tree_bench_double("MC BFS WeightedBinaryTree Benchmark double", 10, 3);

