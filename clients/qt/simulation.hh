/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of HONEI. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef QT_GUARD_SIMULATION_HH
#define QT_GUARD_SIMULATION_HH 1

#include <honei/lbm/tags.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/solver_lbm_grid.hh>

using namespace lbm;
using namespace lbm_lattice_types;


template<typename Tag_,
         typename App_,
         typename Prec_,
         typename ForceScheme_,
         typename ForceType_,
         typename GridType_,
         typename LatticeType_,
         typename BoundaryType_,
         typename Mode_>
class Simulation
{
    private:
        bool _simulation_ready_flag;
        Grid<LatticeType_, Prec_> _grid;
        PackedGridData<LatticeType_, Prec_>  _data;
        PackedGridInfo<LatticeType_> _info;

        SolverLBMGrid<Tag_, App_, Prec_, ForceScheme_, ForceType_, GridType_, LatticeType_, BoundaryType_, Mode_> * _solver;

    public:

        Simulation() :
        _simulation_ready_flag(false)
        {
            ScenarioCollection::get_scenario(0, 50, 50, _grid);
            GridPacker<LatticeType_, BoundaryType_, Prec_>::pack(_grid, _info, _data);

            _solver = new SolverLBMGrid<Tag_, App_, Prec_, ForceScheme_, ForceType_, GridType_, LatticeType_, BoundaryType_, Mode_> (&_info, &_data, _grid.d_x, _grid.d_y, _grid.d_t, _grid.tau);

            _solver->do_preprocessing();

            _simulation_ready_flag = true;
        }

        Simulation(unsigned long scenario_number) :
        _simulation_ready_flag(false)
        {
            ScenarioCollection::get_scenario(scenario_number, 50, 50, _grid);
            GridPacker<LatticeType_, BoundaryType_, Prec_>::pack(_grid, _info, _data);

            _solver = new SolverLBMGrid<Tag_, App_, Prec_, ForceScheme_, ForceType_, GridType_, LatticeType_, BoundaryType_, Mode_> (&_info, &_data, _grid.d_x, _grid.d_y, _grid.d_t, _grid.tau);

            _solver->do_preprocessing();

            _simulation_ready_flag = true;
        }

        Grid<LatticeType_, Prec_> & get_grid()
        {
            return _grid;
        }

        PackedGridData<LatticeType_, Prec_> & get_data()
        {
            return _data;
        }

        PackedGridInfo<LatticeType_> & get_info()
        {
            return _info;
        }

        SolverLBMGrid<Tag_, App_, Prec_, ForceScheme_, ForceType_, GridType_, LatticeType_, BoundaryType_, Mode_> & get_solver()
        {
            return *_solver;
        }

        ~Simulation()
        {
            delete _solver;
            _grid.destroy();
            _data.destroy();
            _info.destroy();
        }
};

#endif
