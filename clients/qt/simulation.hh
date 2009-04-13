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

class DEFAULT;

template<typename Type_>
class Simulation
{
};

template<>
class Simulation<DEFAULT>
{
    private:
        bool _simulation_ready_flag;
        Grid<D2Q9, float> _grid;
        PackedGridData<D2Q9, float>  _data;
        PackedGridInfo<D2Q9> _info;

        SolverLBMGrid<tags::CPU::SSE, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY>* _solver;

    public:

        Simulation() :
        _simulation_ready_flag(false)
        {
            ScenarioCollection::get_scenario(0, 50, 50, _grid);
            GridPacker<D2Q9, NOSLIP, float>::pack(_grid, _info, _data);

            _solver = new SolverLBMGrid<tags::CPU::SSE, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (&_info, &_data, _grid.d_x, _grid.d_y, _grid.d_t, _grid.tau);

            _solver->do_preprocessing();

            _simulation_ready_flag = true;
        }

        Grid<D2Q9, float> & get_grid()
        {
            return _grid;
        }

        PackedGridData<D2Q9, float> & get_data()
        {
            return _data;
        }

        PackedGridInfo<D2Q9> & get_info()
        {
            return _info;
        }

        SolverLBMGrid<tags::CPU::SSE, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> & get_solver()
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
