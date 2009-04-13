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

#ifndef QT_GUARD_SIMULATION_CONTROLLER_HH
#define QT_GUARD_SIMULATION_CONTROLLER_HH 1

#include <simulation.hh>
#include <honei/la/sum.hh>

using namespace lbm;
using namespace lbm_lattice_types;

class SimulationController
{
    private:
        Simulation<DEFAULT> _default_simulation;

    public:
        void do_timestep()
        {
            _default_simulation.get_solver().solve();
            GridPacker<D2Q9, NOSLIP, float>::unpack(_default_simulation.get_grid(), _default_simulation.get_info(), _default_simulation.get_data());
        }

        DenseMatrix<float> & get_h()
        {
            return *_default_simulation.get_grid().h;
        }

        DenseMatrix<float> & get_hb()
        {
            return Sum<tags::CPU::SSE>::value(*_default_simulation.get_grid().h, *_default_simulation.get_grid().b);
        }

        DenseMatrix<float> & get_b()
        {
            return *_default_simulation.get_grid().b;
        }
};

#endif
