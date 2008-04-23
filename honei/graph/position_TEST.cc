/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
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


#include <honei/graph/position.hh>
#include <unittest/unittest.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/util/tags.hh>

#include <string>

using namespace honei;
using  namespace tests;

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
        for (typename MutableMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
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
        for (typename MutableMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
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
        for (typename MutableMatrix<DataType_>::ElementIterator e(Position.begin_elements()),
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
                if ((e.column() == e.row() * 2 +1) || (e.column() == e.row() * 2 +2)) 
                {
                    *e = 1;
                    Edge_Weights[e.column()][e.row()] = 1;
                }
            }
    }
};

template <typename Tag_, typename DataType_>
class WeightedFruchtermanReingoldPositionsQuickTest :
    public QuickTest
{
    public:
        WeightedFruchtermanReingoldPositionsQuickTest(const std::string & type) :
            QuickTest("weighted_fruchterman_reingold_positions_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 1, 
                                8, 1};

            DataType_ node_weights[] =  {2, 2};

            DataType_ edge_weights[] =  {0, 1,
                                         1, 0};

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(2,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++]; 
            }

            i = 0;
            std::tr1::shared_ptr<DenseVector<DataType_> > pNode_Weights(new DenseVector<DataType_>(2));
            for (typename Vector<DataType_>::ElementIterator e(pNode_Weights->begin_elements()),
                    e_end(pNode_Weights->end_elements()); e != e_end ; ++e)
            {
                *e = node_weights[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pEdge_Weights(new SparseMatrix<DataType_>(2,2));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEdge_Weights->begin_elements()),
                    e_end(pEdge_Weights->end_elements()); e != e_end ; ++e)
            {
                if (edge_weights[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = edge_weights[i];
                }
                i++;
            }
            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, *pNode_Weights, *pEdge_Weights);
            position.step_width_factors(1.0, 0.5);

            // update the positions 
            position.update(0.00001,100);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 0.001);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 0.001);
        }
};

WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU, float> weighted_fruchterman_reingold_positions_quick_test_float("float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU, double> weighted_fruchterman_reingold_positions_quick_test_double("double");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore, float> mc_weighted_fruchterman_reingold_positions_quick_test_float("MC float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore, double> mc_weighted_fruchterman_reingold_positions_quick_test_double("MC double");
#ifdef HONEI_SSE
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_fruchterman_reingold_positions_quick_test_float("SSE float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_fruchterman_reingold_positions_quick_test_double("SSE double");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_weighted_fruchterman_reingold_positions_quick_test_float("MC SSE float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_weighted_fruchterman_reingold_positions_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
WeightedFruchtermanReingoldPositionsQuickTest<tags::Cell, float> sse_weighted_fruchterman_reingold_positions_quick_test_float("sse float");
#endif

template <typename Tag_, typename DataType_>
class WeightedKamadaKawaiPositionsQuickTest :
    public QuickTest
{
    public:
        WeightedKamadaKawaiPositionsQuickTest(const std::string & type) :
            QuickTest("weighted_kamada_kawai_positions_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 1,
                                8, 1};

            DataType_ node_weights[] =  {2, 2};

            DataType_ edge_weights[] =  {0, 1,
                                         1, 0};

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(2,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }

            i = 0;
            std::tr1::shared_ptr<DenseVector<DataType_> > pNode_Weights(new DenseVector<DataType_>(2));
            for (typename Vector<DataType_>::ElementIterator e(pNode_Weights->begin_elements()),
                    e_end(pNode_Weights->end_elements()); e != e_end ; ++e)
            {
                *e = node_weights[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pEdge_Weights(new SparseMatrix<DataType_>(2,2));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pEdge_Weights->begin_elements()),
                    e_end(pEdge_Weights->end_elements()); e != e_end ; ++e)
            {
                if (edge_weights[i] > std::numeric_limits<DataType_>::epsilon())
                {
                    *e = edge_weights[i];
                }
                i++;
            }
            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, *pNode_Weights, *pEdge_Weights);

            // update the positions 
            position.update(0.00001,50);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 0.001);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 0.001);
        }
};

WeightedKamadaKawaiPositionsQuickTest<tags::CPU, float> weighted_kamada_kawai_positions_quick_test_float("float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU, double> weighted_kamada_kawai_positions_quick_test_double("double");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::MultiCore, float> mc_weighted_kamada_kawai_positions_quick_test_float("MC float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::MultiCore, double> mc_weighted_kamada_kawai_positions_quick_test_double("MC double");
#ifdef HONEI_SSE
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_kamada_kawai_positions_quick_test_float("SSE float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_kamada_kawai_positions_quick_test_double("SSE double");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_weighted_kamada_kawai_positions_quick_test_float("MC SSE float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_weighted_kamada_kawai_positions_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
WeightedKamadaKawaiPositionsQuickTest<tags::Cell, float> cell_weighted_kamada_kawai_positions_quick_test_float("cell float");
#endif

template <typename Tag_, typename DataType_, typename ST_>
class WeightedFruchtermanReingoldPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        WeightedFruchtermanReingoldPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("weighted_fruchterman_reingold_positions_test<" + type + ">")
        {
            register_tag(Tag_::name);
            if (ST_::get_name() == "Clique") _nodecount = nodecount;
            if (ST_::get_name() == "SquareGrid") _nodecount = nodecount*nodecount;
            if (ST_::get_name() == "BinaryTree") _nodecount = (1 << (nodecount + 1)) -1;
        }

        virtual void run() const
        {
            // Creatoing test scenario
            DenseMatrix<DataType_> Coordinates(_nodecount, 2, DataType_(0));
            DenseVector<DataType_> Node_Weights(_nodecount, DataType_(0));
            SparseMatrix<DataType_> Edge_Weights(_nodecount, _nodecount);
            Scenario<DataType_, ST_>::create(Coordinates, Node_Weights, Edge_Weights);

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(Coordinates, Node_Weights, Edge_Weights);

            // update the positions
            position.update(0.00001, _nodecount * 5);
            std::cout << "number of nodes of WFR  "<< _nodecount << std::endl;
            std::cout << "max_node_force of WFR  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of WFR  "<< position.number_of_iterations() << std::endl;
        }
};

template <typename Tag_, typename DataType_, typename ST_>
class WeightedKamadaKawaiPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        WeightedKamadaKawaiPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("weighted_kamada_kawai_positions_test<" + type + ">")
        {
            register_tag(Tag_::name);
            if (ST_::get_name() == "Clique") _nodecount = nodecount;
            if (ST_::get_name() == "SquareGrid") _nodecount = nodecount*nodecount;
            if (ST_::get_name() == "BinaryTree") _nodecount = (1 << (nodecount + 1)) -1;
        }

        virtual void run() const
        {
            // Creatoing test scenario
            DenseMatrix<DataType_> Coordinates(_nodecount, 2, DataType_(0));
            DenseVector<DataType_> Node_Weights(_nodecount, DataType_(0));
            SparseMatrix<DataType_> Edge_Weights(_nodecount, _nodecount);
            Scenario<DataType_, ST_>::create(Coordinates, Node_Weights, Edge_Weights);

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(Coordinates, Node_Weights, Edge_Weights);

            // update the positions
            //position.statistic(1);
            position.update(0.00001,_nodecount * 10);
            //position.get_statistic("ausgabe.txt");
            std::cout << "number of nodes of WKK  "<< _nodecount << std::endl;
            std::cout << "max_node_force of WKK  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of WKK  "<< position.number_of_iterations() << std::endl;
        }
};

#define POSITIONTEST Scenarios::SquareGrid // possible scenarios are: Clique, SquareGrid, BinaryTree
#define POSITIONTESTSIZE 12 //POSITIONTESTSIZE = numbers of nodes (Clique), POSITIONTESTSIZE = numbers of nodes in a line (SquareGrid), POSITIONTESTSIZE = depth (BinaryTree)

WeightedFruchtermanReingoldPositionsTest<tags::CPU, float, POSITIONTEST> weighted_fruchterman_reingold_positions_test_float("float", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU, double, POSITIONTEST> weighted_fruchterman_reingold_positions_test_double("double", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::MultiCore, float, POSITIONTEST> mc_weighted_fruchterman_reingold_positions_test_float("MC float", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::MultiCore, double, POSITIONTEST> mc_weighted_fruchterman_reingold_positions_test_double("MC double", POSITIONTESTSIZE);
#ifdef HONEI_SSE
WeightedFruchtermanReingoldPositionsTest<tags::CPU::SSE, float, POSITIONTEST> sse_weighted_fruchterman_reingold_positions_test_float("SSE float", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::SSE, double, POSITIONTEST> sse_weighted_fruchterman_reingold_positions_test_double("SSE double", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::MultiCore::SSE, float, POSITIONTEST> mc_sse_weighted_fruchterman_reingold_positions_test_float("MC SSE float", POSITIONTESTSIZE);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::MultiCore::SSE, double, POSITIONTEST> mc_sse_weighted_fruchterman_reingold_positions_test_double("MC SSE double", POSITIONTESTSIZE);
#endif
#ifdef HONEI_CELL
WeightedFruchtermanReingoldPositionsTest<tags::Cell, float, POSITIONTEST> cell_weighted_fruchterman_reingold_positions_test_float("cell float", POSITIONTESTSIZE);
#endif

WeightedKamadaKawaiPositionsTest<tags::CPU, float, POSITIONTEST> weighted_kamada_kawai_positions_test_float("float", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU, double, POSITIONTEST> weighted_kamada_kawai_positions_test_double("double", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU::MultiCore, float, POSITIONTEST> mc_weighted_kamada_kawai_positions_test_float("MC float", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU::MultiCore, double, POSITIONTEST> mc_weighted_kamada_kawai_positions_test_double("MC double", POSITIONTESTSIZE);
#ifdef HONEI_SSE
WeightedKamadaKawaiPositionsTest<tags::CPU::SSE, float, POSITIONTEST> sse_weighted_kamada_kawai_positions_test_float("SSE float", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU::SSE, double, POSITIONTEST> sse_weighted_kamada_kawai_positions_test_double("SSE double", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU::MultiCore::SSE, float, POSITIONTEST> mc_sse_weighted_kamada_kawai_positions_test_float("MC SSE float", POSITIONTESTSIZE);
WeightedKamadaKawaiPositionsTest<tags::CPU::MultiCore::SSE, double, POSITIONTEST> mc_sse_weighted_kamada_kawai_positions_test_double("MC SSE double", POSITIONTESTSIZE);
#endif
#ifdef HONEI_CELL
WeightedKamadaKawaiPositionsTest<tags::Cell, float, POSITIONTEST> cell_weighted_kamada_kawai_positions_test_float("cell float", POSITIONTESTSIZE);
#endif
