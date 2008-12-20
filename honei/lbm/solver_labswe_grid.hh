/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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


#ifndef LBM_GUARD_SOLVER_LABSWE_GRID_HH
#define LBM_GUARD_SOLVER_LABSWE_GRID_HH 1

/**
 * \file
 * Implementation of a SWE solver using LBM and PackedGrid.
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
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sum.hh>
#include <honei/la/scale.hh>
#include <honei/la/element_product.hh>
#include <honei/la/element_inverse.hh>
#include <honei/lbm/collide_stream_grid.hh>
#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/lbm/force_grid.hh>
#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/lbm/extraction_grid.hh>
#include <cmath>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_partitioner.hh>
#include <honei/lbm/grid_packer.hh>

#include <iostream>

using namespace honei::lbm;
using namespace honei::lbm::lbm_boundary_types;

namespace honei
{

    template<typename Tag_,
        typename ResPrec_,
        typename SourceType_,
        typename SourceSheme_,
        typename GridType_,
        typename LatticeType_,
        typename BoundaryType_>
            class SolverLABSWEGrid
            {
            };

                template<typename Tag_, typename ResPrec_>
                    class SolverLABSWEGrid<Tag_, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>
                    {
                        private:
                            /** Global variables.
                             *
                 **/


                ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t;

                unsigned long _time;

                PackedGridData<D2Q9, ResPrec_> * _data;
                PackedGridInfo<D2Q9> * _info;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi;

           public:
                SolverLABSWEGrid(PackedGridData<D2Q9, ResPrec_> * data, PackedGridInfo<D2Q9> * info, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt) :
                    _relaxation_time(ResPrec_(1.5)),
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _time(0),
                    _data(data),
                    _info(info),
                    _n_alpha(ResPrec_(6.)),
                    _gravity(9.80665),
                    _pi(3.14159265)
                    {
                        CONTEXT("When creating LABSWE solver:");
                        _e = _delta_x / _delta_t;

                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                        {
                            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, ResPrec_>::cuda_pack(*_info, *_data);
                        }
                    }

                ~SolverLABSWEGrid()
                {
                    CONTEXT("When destroying LABSWE solver.");
                }

                void do_preprocessing()
                {
                    CONTEXT("When performing LABSWE preprocessing.");

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
                    EquilibriumDistributionGrid<Tag_, lbm_applications::LABSWE>::
                       value(_gravity, _e, *_info, *_data);

                    *_data->f_0 = _data->f_eq_0->copy();
                    *_data->f_1 = _data->f_eq_1->copy();
                    *_data->f_2 = _data->f_eq_2->copy();
                    *_data->f_3 = _data->f_eq_3->copy();
                    *_data->f_4 = _data->f_eq_4->copy();
                    *_data->f_5 = _data->f_eq_5->copy();
                    *_data->f_6 = _data->f_eq_6->copy();
                    *_data->f_7 = _data->f_eq_7->copy();
                    *_data->f_8 = _data->f_eq_8->copy();

                    CollideStreamGrid<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data,
                              _relaxation_time);

#ifdef SOLVER_VERBOSE
                    std::cout << "h after preprocessing:" << std::endl;
                    std::cout << *_data->h << std::endl;
#endif
                }


                /** Capsule for the solution: Single step time marching.
                 *
                 **/
                void solve()
                {
                    ForceGrid<Tag_, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>::value(*_data, *_info, ResPrec_(9.81), _delta_x, _delta_y, _delta_t );

                    ///Boundary correction:
                    UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::
                        value(*_data, *_info);

                    //extract velocities out of h from previous timestep:
                    Extraction<Tag_, lbm_applications::LABSWE>::value(*_info, *_data);

                    ++_time;

                    EquilibriumDistributionGrid<Tag_, lbm_applications::LABSWE>::
                        value(_gravity, _e, *_info, *_data);

                    CollideStreamGrid<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                              *_data,
                              _relaxation_time);
                };
        };

    template<typename ResPrec_>
    class SolverLABSWEGrid<tags::CPU::MultiCore, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>
    {
        private:
            unsigned long _parts;
            PackedGridInfo<D2Q9> * _info;
            PackedGridData<D2Q9, ResPrec_> * _data;
            std::vector<PackedGridInfo<D2Q9> > _info_list;
            std::vector<PackedGridData<D2Q9, ResPrec_> > _data_list;
            std::vector<PackedGridFringe<D2Q9> > _fringe_list;
            std::vector<SolverLABSWEGrid<tags::CPU, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> *> _solver_list;

        public:
            SolverLABSWEGrid(PackedGridData<D2Q9, ResPrec_> * data, PackedGridInfo<D2Q9> * info, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt):
                _parts(4), /// \todo use Configuration
                _info(info),
                _data(data)
        {
            CONTEXT("When creating LABSWE solver:");
            GridPartitioner<D2Q9, ResPrec_>::decompose(_parts, *_info, *_data, _info_list, _data_list, _fringe_list);

            for(unsigned long i(0) ; i < _parts ; ++i)
            {
                _solver_list.push_back(new SolverLABSWEGrid<tags::CPU, ResPrec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>(&_data_list[i], &_info_list[i], 1., 1., 1.));
            }
        }

            ~SolverLABSWEGrid()
            {
                CONTEXT("When destroying LABSWE solver.");
            }

            void do_preprocessing()
            {
                CONTEXT("When performing LABSWE preprocessing.");
                for (unsigned long i(0) ; i < _parts ; ++i)
                {
                    _solver_list.at(i)->do_preprocessing();
                }
                GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
            }

            void solve()
            {
                for (unsigned long i(0) ; i < _parts ; ++i)
                {
                    _solver_list.at(i)->solve();
                }
                GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
                /// \todo remove compose - it is only necessary if one must read the data
                GridPartitioner<D2Q9, ResPrec_>::compose(*_info, *_data, _info_list, _data_list);
            }
    };
}
#endif
