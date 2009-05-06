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

enum solver_type
{
#ifdef HONEI_SSE
    sse_full_dry,
    sse_full_wet,
    sse_wet,
    sse,
#endif
#ifdef HONEI_CUDA
    cuda_full_wet,
    cuda_full_dry,
    cuda_wet,
    cuda,
#endif
    cpu_full_wet,
    cpu_full_dry,
    cpu_wet,
    cpu
};


template <typename Prec_>
class SimulationController
{
    private:
        solver_type _current_solver;

        ///Provide pointers for supported kinds of simulations:
        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cpu_full_dry_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _sse_full_dry_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cuda_full_dry_simulation;
#endif

        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cpu_full_wet_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _sse_full_wet_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cuda_full_wet_simulation;
#endif

        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cpu_wet_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _sse_wet_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cuda_wet_simulation;
#endif

        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cpu_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _sse_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cuda_simulation;
#endif

        bool _noslip_flag;
        bool _d2q9_flag;

    public:
        void do_timestep()
        {
            _sse_full_dry_simulation->get_solver().solve();

            if(_noslip_flag && _d2q9_flag)
                GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_sse_full_dry_simulation->get_grid(), _sse_full_dry_simulation->get_info(), _sse_full_dry_simulation->get_data());
        }

        DenseMatrix<Prec_> & get_h()
        {
            return *_sse_full_dry_simulation->get_grid().h;
        }

        DenseMatrix<Prec_> & get_hb()
        {
#ifdef HONEI_SSE
            return Sum<tags::CPU::SSE>::value(*_sse_full_dry_simulation->get_grid().h, *_sse_full_dry_simulation->get_grid().b);
#else
            return Sum<tags::CPU>::value(*_sse_full_dry_simulation->get_grid().h, *_sse_full_dry_simulation->get_grid().b);
#endif
        }

        DenseMatrix<Prec_> & get_b()
        {
            return *_sse_full_dry_simulation->get_grid().b;
        }

        SimulationController() :
            _noslip_flag(true),
            _d2q9_flag(true),
            _current_solver(sse_full_dry)
        {
            _sse_full_dry_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> ();
        }

        void reload_simulation()
        {
            switch (_current_solver)
            {
                case sse_full_dry:
                    {
                        delete _sse_full_dry_simulation;
                        _sse_full_dry_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> ();
                    }
            }
        }
};

#endif
