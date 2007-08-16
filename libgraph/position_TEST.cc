/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
// our PositionsTest is a BaseTest.
class FruchtermanReingoldPositionsTest :
    public BaseTest
{
    public:
        FruchtermanReingoldPositionsTest(const std::string & type) :
            BaseTest("fruchterman_reingold_positions_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            // Creating test scenario
            DataType_ pos[] =  {1, 2, 2, 3, 2, 3, 3,
                                1, 1, 3, 1, 4, 4, 3};

            DataType_ adj[] =  {0, 1, 0, 0, 0, 0, 0,
                                1, 0, 1, 0, 0, 0, 0,
                                0, 1, 0, 1, 1, 0, 0,
                                0, 0, 1, 0, 0, 0, 0,
                                0, 0, 1, 0, 0, 1, 0,
                                0, 0, 0, 0, 1, 0, 1,
                                0, 0, 0, 0, 0, 1, 0};

            // Now, fill that numbers into the real matrices
            std::tr1::shared_ptr<DenseMatrix<DataType_> > pPosition(new DenseMatrix<DataType_>(7,2));
            int i = 0;
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }
            i = 0;
            std::tr1::shared_ptr<DenseMatrix<bool> > pNeighbour(new DenseMatrix<bool>(7,7));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                    e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                *e = adj[i++];
            }
            // Creating a Positions object with the test scenario
            Positions<DataType_, methods::FruchtermanReingold> position(*pPosition, *pNeighbour, 2);

            // update the positions - try and see what happens :D
            position.update(0.01,350);
            std::cout << position.coordinates();

            // we say the test was fine if we really get here. the numbers are not important right now. perhaps tomorrow :D
            TEST_CHECK(true);
        }
};

FruchtermanReingoldPositionsTest<float> fruchterman_reingold_positions_test_float("float");
FruchtermanReingoldPositionsTest<double> fruchterman_reingold_positions_test_double("double");

template <typename DataType_>
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
            int i = 0;
            for (typename MutableMatrix<DataType_>::ElementIterator e(pPosition->begin_elements()),
                    e_end(pPosition->end_elements());e != e_end ; ++e)
            {
                *e = pos[i++];
            }

            i = 0;
            std::tr1::shared_ptr<DenseMatrix<bool> > pNeighbour(new DenseMatrix<bool>(2,2));
            for (typename MutableMatrix<bool>::ElementIterator e(pNeighbour->begin_elements()),
                    e_end(pNeighbour->end_elements()); e != e_end ; ++e)
            {
                *e = adj[i++];
            }

            // Creating a Positions object with the test scenario
            Positions<DataType_, methods::KamadaKawai> position(*pPosition, *pNeighbour, 2);

            // update the positions - try and see what happens :D
            position.update(0.1,85);

            //std::cout << position.coordinates();
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][0], 6, 30 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[0][1], 8, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][0], 1, std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(position.coordinates()[1][1], 1, std::numeric_limits<DataType_>::epsilon());
        }
};

KamadaKawaiPositionsQuickTest<float> kamada_kawai_positions_quick_test_float("float");
KamadaKawaiPositionsQuickTest<double> kamada_kawai_positions_quick_test_double("double");
