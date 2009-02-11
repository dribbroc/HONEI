/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LBM_GUARD_SCENARIO_COLLECTION_HH
#define LBM_GUARD_SCENARIO_COLLECTION_HH 1

/**
 * \file
 * Implementation of a central scenario database for use in
 * test/benchmarks/clients within the HONEI LBM framework.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/grid.hh>
#include <honei/swe/volume.hh>

using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

class ScenarioCollection
{
    private:
        static const unsigned long _scenario_count = 1;

    public:
        template <typename GridType_, typename DataType_>
        static void get_scenario(unsigned long scen_id, unsigned long grid_height, unsigned long grid_width, Grid<GridType_, DataType_> & target_grid)
        {
            switch (scen_id)
            {
                case 0:
                    {
                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.06), 25, 25);
                        c1.value();
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                        target_grid.description = "Circular dam break, uncritical.\n";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("Discretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                    }
                    break;
            }
            return;
        }

        static unsigned long get_scenario_count()
        {
            return _scenario_count;
        }
};


#endif
