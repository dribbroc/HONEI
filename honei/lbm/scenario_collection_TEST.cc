/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/lbm/scenario_collection.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm::lbm_lattice_types;


template <typename GridType_, typename DataType_>
class ScenarioCollectionTest :
    public BaseTest
{
    public:
        ScenarioCollectionTest(const std::string & type) :
            BaseTest("scenario_collection_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long scen_count(ScenarioCollection::get_scenario_count());
            std::cout << "Scenario count: " << scen_count << std::endl;
            Grid<GridType_, DataType_> grid;

            ScenarioCollection::get_scenario(0, 10, 10, grid);

            std::cout << grid.description << std::endl;
            std::cout << grid.long_description << std::endl;

            std::cout << (*grid.h) << std::endl;

            TEST_CHECK(true);
        }
};
ScenarioCollectionTest<D2Q9, float> sc_test_float("D2Q9, float");
