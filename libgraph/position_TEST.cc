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


#include <libgraph/position.hh>
#include <unittest/unittest.hh>
#include <libla/dense_matrix.hh>
#include <libutil/tags.hh>

#include <string>

using namespace honei;
using  namespace tests;

template <typename Tag_, typename DataType_>
class FruchtermanReingoldPositionsQuickTest :
    public QuickTest
{
    public:
        FruchtermanReingoldPositionsQuickTest(const std::string & type) :
            QuickTest("fruchterman_reingold_positions_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 8, 
                                1, 1};
            DataType_ adj[] = { 0, 1, 
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
            std::tr1::shared_ptr<SparseMatrix<bool> > pNeighbour(new SparseMatrix<bool>(2,2));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            // Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::FruchtermanReingold> position(*pPosition, *pNeighbour, -5);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::FruchtermanReingold> position(*pPosition, *pNeighbour, 2);

            // update the positions 
            position.update(0.01,25);

            //std::cout << "positions of FR  "<< position.coordinates() << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][1] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[1][0], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

FruchtermanReingoldPositionsQuickTest<tags::CPU, float> fruchterman_reingold_positions_quick_test_float("float");
FruchtermanReingoldPositionsQuickTest<tags::CPU, double> fruchterman_reingold_positions_quick_test_double("double");

FruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, float> sse_fruchterman_reingold_positions_quick_test_float("SSE float");
FruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, double> sse_fruchterman_reingold_positions_quick_test_double("SSE double");


template <typename Tag_, typename DataType_>
class KamadaKawaiPositionsQuickTest :
    public QuickTest
{
    public:
        KamadaKawaiPositionsQuickTest(const std::string & type) :
            QuickTest("kamada_kawai_positions_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creatoing test scenario
            DataType_ pos[] = {-2, 8,
                                1, 1 };
            DataType_ adj[] = { 0, 1,
                                1, 0 };

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(2,2));
            int i(0);
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }

            i = 0;
            std::tr1::shared_ptr<SparseMatrix<bool> > pNeighbour(new SparseMatrix<bool>(2,2));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            // Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::FruchtermanReingold> position(*pPosition, *pNeighbour, -5);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }

            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::KamadaKawai> position(*pPosition, *pNeighbour, 2);

            // update the positions 
            position.update(0.01,50);
            // std::cout << "positions of KK  "<< position.coordinates() << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][1] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[1][0], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

KamadaKawaiPositionsQuickTest<tags::CPU, float> kamada_kawai_positions_quick_test_float("float");
KamadaKawaiPositionsQuickTest<tags::CPU, double> kamada_kawai_positions_quick_test_double("double");

//KamadaKawaiPositionsQuickTest<tags::CPU::SSE, float> sse_kamada_kawai_positions_quick_test_float("SSE float");
//KamadaKawaiPositionsQuickTest<tags::CPU::SSE, double> sse_kamada_kawai_positions_quick_test_double("SSE double");

template <typename Tag_, typename DataType_>
class weightedFruchtermanReingoldPositionsQuickTest :
    public QuickTest
{
    public:
        weightedFruchtermanReingoldPositionsQuickTest(const std::string & type) :
            QuickTest("weighted_fruchterman_reingold_positions_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 8, 
                                1, 1};

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
                *e = edge_weights[i++];
            }
            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, *pNode_Weights, *pEdge_Weights);

            // update the positions 
            position.update(0.01,85);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][0], 2, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][1], 4, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0], 1, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1], 1, 85 * std::numeric_limits<DataType_>::epsilon());
        }
};

weightedFruchtermanReingoldPositionsQuickTest<tags::CPU, float> weighted_fruchterman_reingold_positions_quick_test_float("float");
weightedFruchtermanReingoldPositionsQuickTest<tags::CPU, double> weighted_fruchterman_reingold_positions_quick_test_double("double");

weightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_fruchterman_reingold_positions_quick_test_float("SSE float");
weightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_fruchterman_reingold_positions_quick_test_double("SSE double");

template <typename Tag_, typename DataType_>
class weightedKamadaKawaiPositionsQuickTest :
    public QuickTest
{
    public:
        weightedKamadaKawaiPositionsQuickTest(const std::string & type) :
            QuickTest("weighted_kamada_kawai_positions_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 8,
                                1, 1};

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
                *e = edge_weights[i++];
            }
            // Creating a Positions object with the test scenario
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, *pNode_Weights, *pEdge_Weights);

            // update the positions 
            position.update(0.01,85);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][0], 6, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][1], 8, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0], 1, 85 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1], 1, 85 * std::numeric_limits<DataType_>::epsilon());
        }
};

weightedKamadaKawaiPositionsQuickTest<tags::CPU, float> weighted_kamada_kawai_positions_quick_test_float("float");
weightedKamadaKawaiPositionsQuickTest<tags::CPU, double> weighted_kamada_kawai_positions_quick_test_double("double");

weightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_kamada_kawai_positions_quick_test_float("SSE float");
weightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_kamada_kawai_positions_quick_test_double("SSE double");
