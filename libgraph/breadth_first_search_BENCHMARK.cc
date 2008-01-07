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

#include <libgraph/breadth_first_search.hh>


using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class BreadthFirstSearchCliqueBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        BreadthFirstSearchCliqueBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a nodecount-clique
            SparseMatrix<bool>  pAdj(_nodecount, _nodecount);
            for (typename MutableMatrix<bool>::ElementIterator e(pAdj.begin_elements()), e_end(pAdj.end_elements()); e != e_end ; ++e)
            {
                if (e.row() != e.column())
                {
                    *e = true;
                }
            }

            // Creating a distance object with the test scenario
            DenseMatrix<long int> distance(_nodecount, _nodecount, 0);

            // Calculate distance matrix 
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance, pAdj));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class BreadthFirstSearchWeightedCliqueBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        BreadthFirstSearchWeightedCliqueBench(const std::string & id, int nodecount, int count) :
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
            for (typename Vector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(4);
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Calculate distance matrix
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class BreadthFirstSearchBinaryTreeBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        BreadthFirstSearchBinaryTreeBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            // Create a binary tree with deep _nodecount
            SparseMatrix<bool>  pAdj(_nodecount, _nodecount);
            for (unsigned long i(0), column(1); column < _nodecount; i++, column +=2)
            {
                pAdj[i][column] = true;
                pAdj[column][i] = true;
                pAdj[i][column + 1] = true;
                pAdj[column + 1][i] = true;
            }

            // Creating a distance object with the test scenario
            DenseMatrix<long int> distance(_nodecount, _nodecount, 0);

            // Calculate distance matrix
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance, pAdj));
            }

            evaluate();
        }
};

template <typename Tag_, typename DataType_>
class BreadthFirstSearchWeightedBinaryTreeBench :
    public Benchmark
{
    private:
        int _nodecount;
        int _count;
    public:
        BreadthFirstSearchWeightedBinaryTreeBench(const std::string & id, int nodecount, int count) :
            Benchmark(id)
        {
            _nodecount = nodecount;
            _count = count;
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
            for (typename Vector<DataType_>::ElementIterator e(pNW.begin_elements()), e_end(pNW.end_elements()); e != e_end ; ++e)
            {
                *e = DataType_(1 + (exp % 2));
                if (e.index() + 1 == (1 << (exp +1)) -1) exp++;
            }

            // Creating a distance object with the test scenario
            DenseMatrix<DataType_> distance2(_nodecount, _nodecount, DataType_(0));

            // Calculate distance matrix
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(BreadthFirstSearch<Tag_>::value(distance2, pNW, pEW));
            }

            evaluate();
        }
};

BreadthFirstSearchCliqueBench<tags::CPU, float> breadth_first_search_clique_bench_float("BFS Clique Benchmark float", 1000, 3);
BreadthFirstSearchCliqueBench<tags::CPU, double> breadth_first_search_clique_bench_double("BFS Clique Benchmark double", 1000, 3);
BreadthFirstSearchCliqueBench<tags::CPU::MultiCore, float> mc_breadth_first_search_clique_bench_float("MC BFS Clique Benchmark float", 1000, 3);
BreadthFirstSearchCliqueBench<tags::CPU::MultiCore, double> mc_breadth_first_search_clique_bench_double("MC BFS Clique Benchmark double", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU, float> breadth_first_search_weighted_clique_bench_float("BFS WeightedClique Benchmark float", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU, double> breadth_first_search_weighted_clique_bench_double("BFS WeightedClique Benchmark double", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU::MultiCore, float> mc_breadth_first_search_weighted_clique_bench_float("MC BFS WeightedClique Benchmark float", 1000, 3);
BreadthFirstSearchWeightedCliqueBench<tags::CPU::MultiCore, double> mc_breadth_first_search_weighted_clique_bench_double("MC BFS WeightedClique Benchmark double", 1000, 3);
BreadthFirstSearchBinaryTreeBench<tags::CPU, float> breadth_first_search_binary_tree_bench_float("BFS BinaryTree Benchmark float", 30, 3);
BreadthFirstSearchBinaryTreeBench<tags::CPU, double> breadth_first_search_binary_tree_bench_double("BFS BinaryTree Benchmark double", 30, 3);
BreadthFirstSearchBinaryTreeBench<tags::CPU::MultiCore, float> mc_breadth_first_search_binary_tree_bench_float("MC BFS BinaryTree Benchmark float", 30, 3);
BreadthFirstSearchBinaryTreeBench<tags::CPU::MultiCore, double> mc_breadth_first_search_binary_tree_bench_double("MC BFS BinaryTree Benchmark double", 30, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU, float> breadth_first_search_weighted_binary_tree_bench_float("BFS WeightedBinaryTree Benchmark float", 30, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU, double> breadth_first_search_weighted_binary_tree_bench_double("BFS WeightedBinaryTree Benchmark double", 30, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU::MultiCore, float> mc_breadth_first_search_weighted_binary_tree_bench_float("MC BFS WeightedBinaryTree Benchmark float", 30, 3);
BreadthFirstSearchWeightedBinaryTreeBench<tags::CPU::MultiCore, double> mc_breadth_first_search_weighted_binary_tree_bench_double("MC BFS WeightedBinaryTree Benchmark double", 30, 3);

