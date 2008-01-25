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


#include <honei/libgraph/position.hh>
#include <unittest/unittest.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libutil/tags.hh>

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
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] = {-2, 1, 
                                8, 1};
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
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pNeighbour(new SparseMatrix<DataType_>(2,2));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            /*/ Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, nodeWeights, *pNeighbour);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }*/

            // Creating a Positions object with the test scenario
            DenseVector<DataType_> nodeWeights(2, DataType_(1)); 
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, nodeWeights, *pNeighbour);

            // update the positions 
            position.update(0.01,25);

            //std::cout << "positions of FR  "<< position.coordinates() << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

FruchtermanReingoldPositionsQuickTest<tags::CPU, float> fruchterman_reingold_positions_quick_test_float("float");
FruchtermanReingoldPositionsQuickTest<tags::CPU, double> fruchterman_reingold_positions_quick_test_double("double");
//FruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore, float> mc_fruchterman_reingold_positions_quick_test_float("MC float");
//FruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore, double> mc_fruchterman_reingold_positions_quick_test_double("MC double");

#ifdef HONEI_SSE
FruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, float> sse_fruchterman_reingold_positions_quick_test_float("SSE float");
FruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, double> sse_fruchterman_reingold_positions_quick_test_double("SSE double");
//FruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_fruchterman_reingold_positions_quick_test_float("MC SSE float");
//FruchtermanReingoldPositionsQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_fruchterman_reingold_positions_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
FruchtermanReingoldPositionsQuickTest<tags::Cell, float> cell_fruchterman_reingold_positions_quick_test_float("Cell float");
#endif


template <typename Tag_, typename DataType_>
class KamadaKawaiPositionsQuickTest :
    public QuickTest
{
    public:
        KamadaKawaiPositionsQuickTest(const std::string & type) :
            QuickTest("kamada_kawai_positions_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            // Creatoing test scenario
            DataType_ pos[] = {-2, 1,
                                8, 1 };
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
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pNeighbour(new SparseMatrix<DataType_>(2,2));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon())
                {
                    *e = adj[i];
                }
                i++;
            }

            /*/ Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, *pNeighbour, -5);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }*/

            // Creating a Positions object with the test scenario
            DenseVector<DataType_> nodeWeights(2, DataType_(1)); 
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, nodeWeights, *pNeighbour);

            // update the positions 
            position.update(0.01,50);
            // std::cout << "positions of KK  "<< position.coordinates() << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

KamadaKawaiPositionsQuickTest<tags::CPU, float> kamada_kawai_positions_quick_test_float("float");
KamadaKawaiPositionsQuickTest<tags::CPU, double> kamada_kawai_positions_quick_test_double("double");
#ifdef HONEI_SSE
KamadaKawaiPositionsQuickTest<tags::CPU::SSE, float> sse_kamada_kawai_positions_quick_test_float("SSE float");
KamadaKawaiPositionsQuickTest<tags::CPU::SSE, double> sse_kamada_kawai_positions_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
KamadaKawaiPositionsQuickTest<tags::Cell, float> cell_kamada_kawai_positions_quick_test_float("cell_float");
#endif

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

            // update the positions 
            position.update(0.01,50);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU, float> weighted_fruchterman_reingold_positions_quick_test_float("float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU, double> weighted_fruchterman_reingold_positions_quick_test_double("double");
#ifdef HONEI_SSE
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_fruchterman_reingold_positions_quick_test_float("SSE float");
WeightedFruchtermanReingoldPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_fruchterman_reingold_positions_quick_test_double("SSE double");
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
            position.update(0.01,50);
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0] - position.coordinates()[0][0], 2, 50 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1] - position.coordinates()[0][1], 0, 50 * std::numeric_limits<DataType_>::epsilon());
        }
};

WeightedKamadaKawaiPositionsQuickTest<tags::CPU, float> weighted_kamada_kawai_positions_quick_test_float("float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU, double> weighted_kamada_kawai_positions_quick_test_double("double");
#ifdef HONEI_SSE
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, float> sse_weighted_kamada_kawai_positions_quick_test_float("SSE float");
WeightedKamadaKawaiPositionsQuickTest<tags::CPU::SSE, double> sse_weighted_kamada_kawai_positions_quick_test_double("SSE double");
#endif
#ifdef HONEI_CELL
WeightedKamadaKawaiPositionsQuickTest<tags::Cell, float> cell_weighted_kamada_kawai_positions_quick_test_float("cell float");
#endif

template <typename Tag_, typename DataType_>
class KamadaKawaiPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        KamadaKawaiPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("kamada_kawai_positions_test<" + type + ">")
        {
            register_tag(Tag_::name);
            _nodecount = nodecount;
        }

        virtual void run() const
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
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pNeighbour(new SparseMatrix<DataType_>(_nodecount,_nodecount));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon()) 
                {
                    *e = adj[i];
                }
                i++;
            }

            /*/ Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, *pNeighbour, -5);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }*/

            // Creating a Positions object with the test scenario
            DenseVector<DataType_> nodeWeights(_nodecount, DataType_(1)); 
            Positions<Tag_, DataType_, methods::WeightedKamadaKawai> position(*pPosition, nodeWeights, *pNeighbour);

            // update the positions 
            position.update(0.00001, _nodecount * 10);
            std::cout << "max_node_force of KK  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of KK  "<< position.number_of_iterations() << std::endl;
        }
};

KamadaKawaiPositionsTest<tags::CPU, float> kamada_kawai_positions_test_float("float", 100);
KamadaKawaiPositionsTest<tags::CPU, double> kamada_kawai_positions_test_double("double", 100);
#ifdef HONEI_SSE
KamadaKawaiPositionsTest<tags::CPU::SSE, float> sse_kamada_kawai_positions_test_float("SSE float", 100);
KamadaKawaiPositionsTest<tags::CPU::SSE, double> sse_kamada_kawai_positions_test_double("SSE double", 100);
#endif
#ifdef HONEI_CELL
KamadaKawaiPositionsTest<tags::Cell, float> cell_kamada_kawai_positions_test_float("cell float", 100);
#endif

template <typename Tag_, typename DataType_>
class FruchtermanReingoldPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        FruchtermanReingoldPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("fruchterman_reingold_positions_test<" + type + ">")
        {
            register_tag(Tag_::name);
            _nodecount = nodecount;
        }

        virtual void run() const
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
            std::tr1::shared_ptr<SparseMatrix<DataType_> > pNeighbour(new SparseMatrix<DataType_>(_nodecount,_nodecount));
            for (typename MutableMatrix<DataType_>::ElementIterator e(pNeighbour->begin_elements()),
                e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                if (adj[i] > std::numeric_limits<DataType_>::epsilon())
                {
                    *e = adj[i];
                }
                i++;
            }

            /*/ Trying to throw a GraphError
            try
            {
                Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, *pNeighbour, -5);
                TEST_CHECK(false);
            }
            catch (GraphError e)
            {
                TEST_CHECK(true);
            }*/

            // Creating a Positions object with the test scenario
            DenseVector<DataType_> nodeWeights(_nodecount, DataType_(1)); 
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, nodeWeights, *pNeighbour);

            // update the positions 
            position.update(0.00001,_nodecount * 5);
            std::cout << "max_node_force of FR  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of FR  "<< position.number_of_iterations() << std::endl;
        }
};

FruchtermanReingoldPositionsTest<tags::CPU, float> fruchterman_reingold_positions_test_float("float", 100);
FruchtermanReingoldPositionsTest<tags::CPU, double> fruchterman_reingold_positions_test_double("double", 100);
#ifdef HONEI_SSE
FruchtermanReingoldPositionsTest<tags::CPU::SSE, float> sse_fruchterman_reingold_positions_test_float("SSE float", 100);
FruchtermanReingoldPositionsTest<tags::CPU::SSE, double> sse_fruchterman_reingold_positions_test_double("SSE double", 100);
#endif
#ifdef HONEI_CELL
FruchtermanReingoldPositionsTest<tags::Cell, float> cell_fruchterman_reingold_positions_test_float("cell float", 100);
#endif

template <typename Tag_, typename DataType_>
class WeightedFruchtermanReingoldPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        WeightedFruchtermanReingoldPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("weighted_fruchterman_reingold_positions_test<" + type + ">")
        {
            _nodecount = nodecount;
            register_tag(Tag_::name);
        }

        virtual void run() const
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
            Positions<Tag_, DataType_, methods::WeightedFruchtermanReingold> position(*pPosition, *pNode_Weights, *pEdge_Weights);

            // update the positions
            position.update(0.00001,_nodecount * 5);
            std::cout << "max_node_force of WFR  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of WFR  "<< position.number_of_iterations() << std::endl;
        }
};

WeightedFruchtermanReingoldPositionsTest<tags::CPU, float> weighted_fruchterman_reingold_positions_test_float("float", 100);
WeightedFruchtermanReingoldPositionsTest<tags::CPU, double> weighted_fruchterman_reingold_positions_test_double("double", 100);
#ifdef HONEI_SSE
WeightedFruchtermanReingoldPositionsTest<tags::CPU::SSE, float> sse_weighted_fruchterman_reingold_positions_test_float("SSE float", 100);
WeightedFruchtermanReingoldPositionsTest<tags::CPU::SSE, double> sse_weighted_fruchterman_reingold_positions_test_double("SSE double", 100);
#endif
#ifdef HONEI_CELL
WeightedFruchtermanReingoldPositionsTest<tags::Cell, float> cell_weighted_fruchterman_reingold_positions_test_float("cell float", 100);
#endif

template <typename Tag_, typename DataType_>
class WeightedKamadaKawaiPositionsTest :
    public BaseTest
{
    private:

        int _nodecount;

    public:
        WeightedKamadaKawaiPositionsTest(const std::string & type, int nodecount = 50) :
            BaseTest("weighted_kamada_kawai_positions_test<" + type + ">")
        {
            _nodecount = nodecount;
            register_tag(Tag_::name);
        }

        virtual void run() const
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
            position.update(0.00001,_nodecount * 10);
            std::cout << "max_node_force of WKK  "<< position.max_node_force() << std::endl;
            std::cout << "number_of_iterations of WKK  "<< position.number_of_iterations() << std::endl;
        }
};
WeightedKamadaKawaiPositionsTest<tags::CPU, float> weighted_kamada_kawai_positions_test_float("float", 100);
WeightedKamadaKawaiPositionsTest<tags::CPU, double> weighted_kamada_kawai_positions_test_double("double", 100);
#ifdef HONEI_SSE
WeightedKamadaKawaiPositionsTest<tags::CPU::SSE, float> sse_weighted_kamada_kawai_positions_test_float("SSE float", 100);
WeightedKamadaKawaiPositionsTest<tags::CPU::SSE, double> sse_weighted_kamada_kawai_positions_test_double("SSE double", 100);
#endif
#ifdef HONEI_CELL
WeightedKamadaKawaiPositionsTest<tags::Cell, float> cell_weighted_kamada_kawai_positions_test_float("cell float", 100);
#endif
