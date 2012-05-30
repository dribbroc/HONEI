/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <apryde@gmx.de>
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


#pragma once
#ifndef LBM_GUARD_SOLVER_LBM_GRID_POLLUTANT_HH
#define LBM_GUARD_SOLVER_LBM_GRID_POLLUTANT_HH 1

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
#include <honei/lbm/collide_stream_grid_pollutant.hh>
#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/lbm/equilibrium_distribution_grid_pollutant.hh>
#include <honei/lbm/force_grid.hh>
#include <honei/lbm/force_grid_pollutant.hh>
#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/lbm/update_velocity_directions_grid_pollutant.hh>
#include <honei/lbm/extraction_grid.hh>
#include <honei/lbm/extraction_grid_pollutant.hh>
#include <cmath>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_partitioner.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/backends/multicore/dispatch_policy.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/backends/cuda/transfer.hh>

#include <iostream>
#include <honei/util/tr1_boost.hh>
#include <vector>


using namespace honei::lbm;
using namespace honei::lbm::lbm_boundary_types;

namespace honei
{
    class SolverLBMGridPollutantBase
    {
        public:
            virtual void do_preprocessing() = 0;
            virtual void do_postprocessing() = 0;
            virtual void solve() = 0;

    };

    template<typename Tag_,
        typename Application_,
        typename ResPrec_,
        typename Force_,
        typename SourceScheme_,
        typename GridType_,
        typename LatticeType_,
        typename BoundaryType_,
        typename LbmMode_>
            class SolverLBMGridPollutant
            {
            };

    template<typename Tag_, typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMGridPollutant<Tag_, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> : public SolverLBMGridBase
        {
            private:
                /** Global variables.
                 *
                 **/


                ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t, _k, _s_0, _dir_value;

                unsigned long _time;

                PackedGridInfo<D2Q9> * _info;
                PackedGridData<D2Q9, ResPrec_> * _data_flow;
                PackedGridData<D2Q9, ResPrec_> * _data_poll;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi, _e_squared;

            public:
                SolverLBMGridPollutant(PackedGridInfo<D2Q9> * info,
                        PackedGridData<D2Q9, ResPrec_> * data_flow,
                        PackedGridData<D2Q9, ResPrec_> * data_poll,
                        ResPrec_ dx,
                        ResPrec_ dy,
                        ResPrec_ dt,
                        ResPrec_ rel_time,
                        ResPrec_ k,
                        ResPrec_ s_0,
                        ResPrec_ dirvalue) :
                    _relaxation_time(rel_time),
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _k(k),
                    _s_0(s_0),
                    _dir_value(dirvalue),
                    _time(0),
                    _info(info),
                    _data_flow(data_flow),
                    _data_poll(data_poll),
                    _n_alpha(ResPrec_(6.)),
                    _gravity(9.80665),
                    _pi(3.14159265)
            {
                CONTEXT("When creating solver for pollutant transport:");
                _e = _delta_x / _delta_t;
                _e_squared = _e * _e;

                if (Tag_::tag_value == tags::tv_gpu_cuda)
                {
                    GridPacker<D2Q9, lbm_boundary_types::NOSLIP, ResPrec_>::cuda_pack(*_info, *_data_flow);
                }
            }

                virtual ~SolverLBMGridPollutant()
                {
                    CONTEXT("When destroying LABSWE solver.");
                }

                void do_preprocessing()
                {
                    CONTEXT("When performing LABSWE preprocessing.");

                    (*_data_poll->distribution_x)[0] = ResPrec_(0.);
                    (*_data_poll->distribution_x)[1] = ResPrec_(1.);
                    (*_data_poll->distribution_x)[2] = ResPrec_(1.);
                    (*_data_poll->distribution_x)[3] = ResPrec_(0.);
                    (*_data_poll->distribution_x)[4] = ResPrec_(-1.);
                    (*_data_poll->distribution_x)[5] = ResPrec_(-1);
                    (*_data_poll->distribution_x)[6] = ResPrec_(-1);
                    (*_data_poll->distribution_x)[7] = ResPrec_(0);
                    (*_data_poll->distribution_x)[8] = ResPrec_(1);
                    (*_data_poll->distribution_y)[0] = ResPrec_(0.);
                    (*_data_poll->distribution_y)[1] = ResPrec_(0.);
                    (*_data_poll->distribution_y)[2] = ResPrec_(1.);
                    (*_data_poll->distribution_y)[3] = ResPrec_(1.);
                    (*_data_poll->distribution_y)[4] = ResPrec_(1.);
                    (*_data_poll->distribution_y)[5] = ResPrec_(0);
                    (*_data_poll->distribution_y)[6] = ResPrec_(-1);
                    (*_data_poll->distribution_y)[7] = ResPrec_(-1);
                    (*_data_poll->distribution_y)[8] = ResPrec_(-1);

                    ///Compute initial equilibrium distribution:
                    EquilibriumDistributionGridPollutant<Tag_, Application_>::
                        value(_e_squared, *_info, *_data_flow, *_data_poll);

                    *_data_poll->f_0 = _data_poll->f_eq_0->copy();
                    *_data_poll->f_1 = _data_poll->f_eq_1->copy();
                    *_data_poll->f_2 = _data_poll->f_eq_2->copy();
                    *_data_poll->f_3 = _data_poll->f_eq_3->copy();
                    *_data_poll->f_4 = _data_poll->f_eq_4->copy();
                    *_data_poll->f_5 = _data_poll->f_eq_5->copy();
                    *_data_poll->f_6 = _data_poll->f_eq_6->copy();
                    *_data_poll->f_7 = _data_poll->f_eq_7->copy();
                    *_data_poll->f_8 = _data_poll->f_eq_8->copy();

                    CollideStreamGridPollutant<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data_flow,
                              *_data_poll,
                              _relaxation_time);

#ifdef SOLVER_VERBOSE
                    std::cout << "h after preprocessing:" << std::endl;
                    std::cout << *_data_flow->h << std::endl;
#endif
                }

                void do_postprocessing()
                {
                }


                /** Capsule for the solution: Single step time marching.
                 *
                 **/
                void solve()
                {
                    ForceGridPollutant<Tag_, Application_>::value(*_info, *_data_flow, *_data_poll, _delta_t, _k, _s_0);

                    ///Boundary correction:
                    //UpdateVelocityDirectionsGridPollutant<Tag_, DIRICHLET_SLIP>::
                      //  value(*_info, *_data_flow, *_data_poll, _e_squared, _dir_value);
                    UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::
                        value(*_info, *_data_poll);

                    //extract velocities out of h from previous timestep:
                    ExtractionGridPollutant<Tag_, LbmMode_>::value(*_info, *_data_flow, *_data_poll, ResPrec_(10e-5));

                    ++_time;

                    EquilibriumDistributionGridPollutant<Tag_, Application_>::
                        value(_e_squared, *_info, *_data_flow, *_data_poll);

                    CollideStreamGridPollutant<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data_flow,
                              *_data_poll,
                              _relaxation_time);
                }
        };
}
#endif
