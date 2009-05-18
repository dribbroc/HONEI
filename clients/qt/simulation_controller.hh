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
#endif
#ifdef HONEI_CUDA
    cuda_full_wet,
    cuda_full_dry,
    cuda_wet,
#endif
    cpu_full_wet,
    cpu_full_dry,
    cpu_wet
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
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> * _cuda_full_dry_simulation;
#endif

        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _cpu_full_wet_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _sse_full_wet_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _cuda_full_wet_simulation;
#endif

        Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _cpu_wet_simulation;

#ifdef HONEI_SSE
        Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _sse_wet_simulation;
#endif

#ifdef HONEI_CUDA
        Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> * _cuda_wet_simulation;
#endif

        bool _noslip_flag;
        bool _d2q9_flag;

    public:
        ///Constructor
        SimulationController() :
            _current_solver(cpu_full_dry), ///Init wit CPU
            _noslip_flag(true),
            _d2q9_flag(true)
    {
        _cpu_full_dry_simulation = new Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> ();

#ifdef HONEI_SSE ///Reset to SSE if possible
        _sse_full_dry_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> ();
        _sse_full_wet_simulation = 0;
        _sse_wet_simulation = 0;

        _current_solver = sse_full_dry;
        delete _cpu_full_dry_simulation;
        _cpu_full_dry_simulation = 0;
#endif
#ifdef HONEI_CUDA
        _cuda_full_dry_simulation = 0;
        _cuda_full_wet_simulation = 0;
        _cuda_wet_simulation = 0;
#endif
        _cpu_full_wet_simulation = 0;
        _cpu_wet_simulation = 0;
    }

        void do_timestep()
        {
            switch(_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        _sse_full_dry_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_sse_full_dry_simulation->get_grid(), _sse_full_dry_simulation->get_info(), _sse_full_dry_simulation->get_data());
                    }
                    break;
                case sse_full_wet:
                    {
                        _sse_full_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_sse_full_wet_simulation->get_grid(), _sse_full_wet_simulation->get_info(), _sse_full_wet_simulation->get_data());
                    }
                    break;
                case sse_wet:
                    {
                        _sse_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_sse_wet_simulation->get_grid(), _sse_wet_simulation->get_info(), _sse_wet_simulation->get_data());
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        _cuda_full_dry_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, float>::unpack(_cuda_full_dry_simulation->get_grid(), _cuda_full_dry_simulation->get_info(), _cuda_full_dry_simulation->get_data());
                    }
                    break;
                case cuda_full_wet:
                    {
                        _cuda_full_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, float>::unpack(_cuda_full_wet_simulation->get_grid(), _cuda_full_wet_simulation->get_info(), _cuda_full_wet_simulation->get_data());
                    }
                    break;
                case cuda_wet:
                    {
                        _cuda_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, float>::unpack(_cuda_wet_simulation->get_grid(), _cuda_wet_simulation->get_info(), _cuda_wet_simulation->get_data());
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        _cpu_full_dry_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_cpu_full_dry_simulation->get_grid(), _cpu_full_dry_simulation->get_info(), _cpu_full_dry_simulation->get_data());
                    }
                    break;
                case cpu_full_wet:
                    {
                        _cpu_full_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_cpu_full_wet_simulation->get_grid(), _cpu_full_wet_simulation->get_info(), _cpu_full_wet_simulation->get_data());
                    }
                    break;
                case cpu_wet:
                    {
                        _cpu_wet_simulation->get_solver().solve();

                        if(_noslip_flag && _d2q9_flag)
                            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_cpu_wet_simulation->get_grid(), _cpu_wet_simulation->get_info(), _cpu_wet_simulation->get_data());
                    }
                    break;

            }
        }

        DenseMatrix<Prec_> & get_h()
        {
            switch(_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        return *_sse_full_dry_simulation->get_grid().h;
                    }
                    break;
                case sse_full_wet:
                    {
                        return *_sse_full_wet_simulation->get_grid().h;
                    }
                    break;
                case sse_wet:
                    {
                        return *_sse_wet_simulation->get_grid().h;
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        return *_cuda_full_dry_simulation->get_grid().h;
                    }
                    break;
                case cuda_full_wet:
                    {
                        return *_cuda_full_wet_simulation->get_grid().h;
                    }
                    break;
                case cuda_wet:
                    {
                        return *_cuda_wet_simulation->get_grid().h;
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        return *_cpu_full_dry_simulation->get_grid().h;
                    }
                    break;
                case cpu_full_wet:
                    {
                        return *_cpu_full_wet_simulation->get_grid().h;
                    }
                    break;
                case cpu_wet:
                    {
                        return *_cpu_wet_simulation->get_grid().h;
                    }
                    break;
            }
        }

        DenseMatrix<Prec_> get_hb()
        {
            switch(_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        DenseMatrix<Prec_> result((*_sse_full_dry_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_sse_full_dry_simulation->get_grid().b));
                        return result;
                    }
                    break;
                case sse_full_wet:
                    {
                        DenseMatrix<Prec_> result((*_sse_full_wet_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_sse_full_wet_simulation->get_grid().b));
                        return result;
                    }
                    break;
                case sse_wet:
                    {
                        DenseMatrix<Prec_> result((*_sse_wet_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_sse_wet_simulation->get_grid().b));
                        return result;
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        DenseMatrix<float> result((*_cuda_full_dry_simulation->get_grid().h).copy());
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());

                        Sum<tags::GPU::CUDA>::value(result, (*_cuda_full_dry_simulation->get_grid().b));

                        convert(converted_result, result);

                        return converted_result;
                    }
                    break;
                case cuda_full_wet:
                    {
                        DenseMatrix<float> result((*_cuda_full_wet_simulation->get_grid().h).copy());
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());

                        Sum<tags::GPU::CUDA>::value(result, (*_cuda_full_wet_simulation->get_grid().b));

                        convert(converted_result, result);

                        return converted_result;
                    }
                    break;
                case cuda_wet:
                    {
                        DenseMatrix<float> result((*_cuda_wet_simulation->get_grid().h).copy());
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());

                        Sum<tags::GPU::CUDA>::value(result, (*_cuda_wet_simulation->get_grid().b));

                        convert(converted_result, result);

                        return converted_result;
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        DenseMatrix<Prec_> result((*_cpu_full_dry_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_cpu_full_dry_simulation->get_grid().b));
                        return result;
                    }
                    break;
                case cpu_full_wet:
                    {
                        DenseMatrix<Prec_> result((*_cpu_full_wet_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_cpu_full_wet_simulation->get_grid().b));
                        return result;
                    }
                    break;
                case cpu_wet:
                    {
                        DenseMatrix<Prec_> result((*_cpu_wet_simulation->get_grid().h).copy());
                        Sum<tags::CPU::SSE>::value(result, (*_cpu_wet_simulation->get_grid().b));
                        return result;
                    }
                    break;
                default:
                    throw InternalError("Undefined solver state.");
            }
        }

        DenseMatrix<Prec_> get_b()
        {
            switch(_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        return *_sse_full_dry_simulation->get_grid().b;
                    }
                    break;
                case sse_full_wet:
                    {
                        return *_sse_full_wet_simulation->get_grid().b;
                    }
                    break;
                case sse_wet:
                    {
                        return *_sse_wet_simulation->get_grid().b;
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        DenseMatrix<float> result(*_cuda_full_dry_simulation->get_grid().b);
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());
                        convert(converted_result, result);
                        return converted_result;
                    }
                    break;
                case cuda_full_wet:
                    {
                        DenseMatrix<float> result(*_cuda_full_wet_simulation->get_grid().b);
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());
                        convert(converted_result, result);
                        return converted_result;
                    }
                    break;
                case cuda_wet:
                    {
                        DenseMatrix<float> result(*_cuda_wet_simulation->get_grid().b);
                        DenseMatrix<Prec_> converted_result(result.rows(), result.columns());
                        convert(converted_result, result);
                        return converted_result;
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        return *_cpu_full_dry_simulation->get_grid().b;
                    }
                    break;
                case cpu_full_wet:
                    {
                        return *_cpu_full_wet_simulation->get_grid().b;
                    }
                    break;
                case cpu_wet:
                    {
                        return *_cpu_wet_simulation->get_grid().b;
                    }
                    break;
                default:
                    throw InternalError("Undefined solver state.");
            }
        }

        void load_simulation(unsigned long id)
        {
            switch (_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        delete _sse_full_dry_simulation;
                        _sse_full_dry_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (id);
                    }
                case sse_full_wet:
                    {
                        delete _sse_full_wet_simulation;
                        _sse_full_wet_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
                case sse_wet:
                    {
                        delete _sse_wet_simulation;
                        _sse_wet_simulation = new Simulation<tags::CPU::SSE, lbm_applications::LABSWE, Prec_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        delete _cuda_full_dry_simulation;
                        _cuda_full_dry_simulation = new Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (id);
                    }
                case cuda_full_wet:
                    {
                        delete _cuda_full_wet_simulation;
                        _cuda_full_wet_simulation = new Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
                case cuda_wet:
                    {
                        delete _cuda_wet_simulation;
                        _cuda_wet_simulation = new Simulation<tags::GPU::CUDA, lbm_applications::LABSWE, float,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
#endif
                case cpu_full_dry:
                    {
                        delete _cpu_full_dry_simulation;
                        _cpu_full_dry_simulation = new Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (id);
                    }
                case cpu_full_wet:
                    {
                        delete _cpu_full_wet_simulation;
                        _cpu_full_wet_simulation = new Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
                case cpu_wet:
                    {
                        delete _cpu_wet_simulation;
                        _cpu_wet_simulation = new Simulation<tags::CPU, lbm_applications::LABSWE, Prec_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> (id);
                    }
            }
        }


        std::string get_backend_info()
        {
            switch (_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        return "SSE";
                    }
                    break;
                case sse_full_wet:
                    {
                        return "SSE";
                    }
                    break;
                case sse_wet:
                    {
                        return "SSE";
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        return "CUDA";
                    }
                    break;
                case cuda_full_wet:
                    {
                        return "CUDA";
                    }
                    break;
                case cuda_wet:
                    {
                        return "CUDA";
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        return "CPU";
                    }
                    break;
                case cpu_full_wet:
                    {
                        return "CPU";
                    }
                    break;
                case cpu_wet:
                    {
                        return "CPU";
                    }
                    break;
                default:
                    throw InternalError("Undefined solver state.");
            }
        }

        std::string get_simulation_info()
        {
            switch (_current_solver)
            {
#ifdef HONEI_SSE
                case sse_full_dry:
                    {
                        return _sse_full_dry_simulation->get_grid().description;
                    }
                    break;
                case sse_full_wet:
                    {
                        return _sse_full_wet_simulation->get_grid().description;
                    }
                    break;
                case sse_wet:
                    {
                        return _sse_wet_simulation->get_grid().description;
                    }
                    break;
#endif
#ifdef HONEI_CUDA
                case cuda_full_dry:
                    {
                        return _cuda_full_dry_simulation->get_grid().description;
                    }
                    break;
                case cuda_full_wet:
                    {
                        return _cuda_full_wet_simulation->get_grid().description;
                    }
                    break;
                case cuda_wet:
                    {
                        return _cuda_wet_simulation->get_grid().description;
                    }
                    break;
#endif
                case cpu_full_dry:
                    {
                        return _cpu_full_dry_simulation->get_grid().description;
                    }
                    break;
                case cpu_full_wet:
                    {
                        return _cpu_full_wet_simulation->get_grid().description;
                    }
                    break;
                case cpu_wet:
                    {
                        return _cpu_wet_simulation->get_grid().description;
                    }
                    break;
                default:
                    throw InternalError("Undefined solver state.");
            }
        }
};

#endif
