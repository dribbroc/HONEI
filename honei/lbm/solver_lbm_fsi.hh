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


#ifndef LBM_GUARD_SOLVER_LBM_GRID_HH
#define LBM_GUARD_SOLVER_LBM_GRID_HH 1

/**
 * \file
 * Implementation of a SWE solver using LBM and PackedGrid with FSI.
 *
 * \ingroup grpliblbm
 **/

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

#ifdef DEBUG
#define SOLVER_VERBOSE
#endif

#include <honei/lbm/tags.hh>
#include <honei/util/tags.hh>
#include <honei/util/configuration.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sum.hh>
#include <honei/la/scale.hh>
#include <honei/la/element_product.hh>
#include <honei/la/element_inverse.hh>
#include <honei/lbm/collide_stream_grid.hh>
#include <honei/lbm/collide_stream_fsi.hh>
#include <honei/lbm/boundary_init_fsi.hh>
#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/lbm/force_grid.hh>
#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/lbm/extraction_grid.hh>
#include <cmath>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_partitioner.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/backends/multicore/dispatch_policy.hh>
#include <honei/backends/multicore/thread_pool.hh>

#include <iostream>
#include <tr1/functional>
#include <vector>

using namespace honei::lbm;
using namespace honei::lbm::lbm_boundary_types;

namespace honei
{

    template<typename Tag_,
        typename Application_,
        typename ResPrec_,
        typename Force_,
        typename SourceScheme_,
        typename GridType_,
        typename LatticeType_,
        typename BoundaryType_,
        typename LbmMode_>
            class SolverLBMFSI
            {
            };

    template<typename Tag_, typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMFSI<Tag_, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>
        {
            private:
                /** Global variables.
                 *
                 **/
                ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t;

                unsigned long _time;

                PackedGridInfo<D2Q9> * _info;
                PackedGridData<D2Q9, ResPrec_> * _data;
                PackedSolidData<D2Q9, ResPrec_> * _solids;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi, _e_squared;

            public:
                SolverLBMFSI(PackedGridInfo<D2Q9> * info,
                             PackedGridData<D2Q9, ResPrec_> * data,
                             PackedSolidData<D2Q9, ResPrec_> * solids,
                             ResPrec_ dx,
                             ResPrec_ dy,
                             ResPrec_ dt,
                             ResPrec_ rel_time) :
                    _relaxation_time(rel_time),
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _time(0),
                    _info(info),
                    _data(data),
                    _solids(solids),
                    _n_alpha(ResPrec_(6.)),
                    _gravity(9.80665),
                    _pi(3.14159265)
            {
                CONTEXT("When creating LBM FSI solver:");
                _e = _delta_x / _delta_t;
                _e_squared = _e * _e;

                /*if (Tag_::tag_value == tags::tv_gpu_cuda)
                {
                    GridPacker<D2Q9, lbm_boundary_types::NOSLIP, ResPrec_>::cuda_pack(*_info, *_data);
                }*/
            }

                ~SolverLBMFSI()
                {
                    CONTEXT("When destroying LBM FSI solver.");
                }

                void do_preprocessing()
                {
                    CONTEXT("When performing LBM FSI preprocessing.");

                    (*_data->distribution_x)[0] = ResPrec_(0.);
                    (*_data->distribution_x)[1] = ResPrec_(_e * cos(ResPrec_(0.)));
                    (*_data->distribution_x)[2] = ResPrec_(sqrt(ResPrec_(2.)) * _e * cos(_pi / ResPrec_(4.)));
                    (*_data->distribution_x)[3] = ResPrec_(_e * cos(_pi / ResPrec_(2.)));
                    (*_data->distribution_x)[4] = ResPrec_(sqrt(ResPrec_(2.)) * _e * cos(ResPrec_(3.) * _pi / ResPrec_(4.)));
                    (*_data->distribution_x)[5] = ResPrec_(_e * cos(_pi));
                    (*_data->distribution_x)[6] = ResPrec_(sqrt(ResPrec_(2.)) * _e * cos(ResPrec_(5.) * _pi / ResPrec_(4.)));
                    (*_data->distribution_x)[7] = ResPrec_(_e * cos(ResPrec_(3.) * _pi / ResPrec_(2.)));
                    (*_data->distribution_x)[8] = ResPrec_(sqrt(ResPrec_(2.)) * _e * cos(ResPrec_(7.) * _pi / ResPrec_(4.)));
                    (*_data->distribution_y)[0] = ResPrec_(0.);
                    (*_data->distribution_y)[1] = ResPrec_(_e * sin(ResPrec_(0.)));
                    (*_data->distribution_y)[2] = ResPrec_(sqrt(ResPrec_(2.)) * _e * sin(_pi / ResPrec_(4.)));
                    (*_data->distribution_y)[3] = ResPrec_(_e * sin(_pi / ResPrec_(2.)));
                    (*_data->distribution_y)[4] = ResPrec_(sqrt(ResPrec_(2.)) * _e * sin(ResPrec_(3.) * _pi / ResPrec_(4.)));
                    (*_data->distribution_y)[5] = ResPrec_(_e * sin(_pi));
                    (*_data->distribution_y)[6] = ResPrec_(sqrt(ResPrec_(2.)) * _e * sin(ResPrec_(5.) * _pi / ResPrec_(4.)));
                    (*_data->distribution_y)[7] = ResPrec_(_e * sin(ResPrec_(3.) * _pi / ResPrec_(2.)));
                    (*_data->distribution_y)[8] = ResPrec_(sqrt(ResPrec_(2.)) * _e * sin(ResPrec_(7.) * _pi / ResPrec_(4.)));

                    ///Compute initial equilibrium distribution:
                    EquilibriumDistributionGrid<Tag_, Application_>::
                        value(_gravity, _e_squared, *_info, *_data);

                    *_data->f_0 = _data->f_eq_0->copy();
                    *_data->f_1 = _data->f_eq_1->copy();
                    *_data->f_2 = _data->f_eq_2->copy();
                    *_data->f_3 = _data->f_eq_3->copy();
                    *_data->f_4 = _data->f_eq_4->copy();
                    *_data->f_5 = _data->f_eq_5->copy();
                    *_data->f_6 = _data->f_eq_6->copy();
                    *_data->f_7 = _data->f_eq_7->copy();
                    *_data->f_8 = _data->f_eq_8->copy();

                    CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data,
                              _relaxation_time);
                    CollideStreamFSI<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info, *_data, *_solids, _delta_x, _delta_y);

#ifdef SOLVER_VERBOSE
                    std::cout << "h after preprocessing:" << std::endl;
                    std::cout << *_data->h << std::endl;
#endif
                }

                void do_postprocessing()
                {
                }


                /** Capsule for the solution: Single step time marching.
                 *
                 **/
                void solve(unsigned long dir)
                {
                    ForceGrid<Tag_, Application_, Force_, SourceScheme_>::value(*_info, *_data, ResPrec_(9.81), _delta_x, _delta_y, _delta_t, ResPrec_(0.01));

                    ///Boundary correction:
                    UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::
                        value(*_info, *_data, _relaxation_time);

                    //extract velocities out of h from previous timestep:
                    ExtractionGrid<Tag_, LbmMode_>::value(*_info, *_data, ResPrec_(10e-5));

                    ++_time;

                    switch(dir)
                    {
                        case 1:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_1>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 2:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_2>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 3:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_3>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 4:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_4>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 5:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_5>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 6:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_6>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 7:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_7>::value(*_info, *_data, *_solids);
                            }
                            break;
                        case 8:
                            {
                                BoundaryInitFSI<Tag_, D2Q9::DIR_8>::value(*_info, *_data, *_solids);
                            }
                            break;
                    }
                    EquilibriumDistributionGrid<Tag_, Application_>::
                        value(_gravity, _e_squared, *_info, *_data);

                    CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data,
                              _relaxation_time);
                    CollideStreamFSI<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info, *_data, *_solids, _delta_x, _delta_y);
                }

        };
}
#endif
