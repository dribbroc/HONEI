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


#pragma once
#ifndef LBM_GUARD_SOLVER_LBM_GRID_HH
#define LBM_GUARD_SOLVER_LBM_GRID_HH 1

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
#include <honei/util/configuration.hh>
#include <honei/util/benchmark_info.hh>
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
    class SolverLBMGridBase
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
            class SolverLBMGrid
            {
            };

    template<typename Tag_, typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMGrid<Tag_, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> : public SolverLBMGridBase
        {
            private:
                /** Global variables.
                 *
                 **/


                ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t;

                unsigned long _time;

                PackedGridInfo<D2Q9> * _info;
                PackedGridData<D2Q9, ResPrec_> * _data;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi, _e_squared;

            public:
                SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time) :
                    _relaxation_time(rel_time),
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _time(0),
                    _info(info),
                    _data(data),
                    _n_alpha(ResPrec_(6.)),
                    _gravity(9.80665),
                    _pi(3.14159265)
            {
                CONTEXT("When creating LABSWE solver:");
                _e = _delta_x / _delta_t;
                _e_squared = _e * _e;

                if (Tag_::tag_value == tags::tv_gpu_cuda)
                {
                    GridPacker<D2Q9, lbm_boundary_types::NOSLIP, ResPrec_>::cuda_pack(*_info, *_data);
                }
            }

                virtual ~SolverLBMGrid()
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
                void solve()
                {
                    ForceGrid<Tag_, Application_, Force_, SourceScheme_>::value(*_info, *_data, ResPrec_(9.81), _delta_x, _delta_y, _delta_t, ResPrec_(0.01));

                    ///Boundary correction:
                    UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::
                        value(*_info, *_data);

                    //extract velocities out of h from previous timestep:
                    ExtractionGrid<Tag_, LbmMode_>::value(*_info, *_data, ResPrec_(10e-5));

                    ++_time;

                    EquilibriumDistributionGrid<Tag_, Application_>::
                        value(_gravity, _e_squared, *_info, *_data);

                    CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                        value(*_info,
                                *_data,
                                _relaxation_time);
                }

                static LBMBenchmarkInfo get_benchmark_info(Grid<D2Q9, ResPrec_> * grid, PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data)
                {
                    LBMBenchmarkInfo result;
                    BenchmarkInfo eq_dist(EquilibriumDistributionGrid<Tag_, Application_>::get_benchmark_info(info, data));
                    result += eq_dist;
                    BenchmarkInfo col_stream(CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::get_benchmark_info(info, data));
                    result += col_stream;
                    BenchmarkInfo force(ForceGrid<Tag_, Application_, Force_, SourceScheme_>::get_benchmark_info(info, data));
                    result += force;
                    BenchmarkInfo extraction(ExtractionGrid<Tag_, LbmMode_>::get_benchmark_info(info, data));
                    result += extraction;

                    result.size.push_back(grid->h->rows());
                    result.size.push_back(grid->h->columns());
                    result.lups = grid->h->rows() * grid->h->columns();
                    result.flups = data->h->size();
                    return result;
                }
        };

    namespace mc
    {
        template<typename Tag_,
            typename Application_,
            typename ResPrec_,
            typename Force__,
            typename SourceScheme_,
            typename GridType_,
            typename LatticeType_,
            typename BoundaryType_,
            typename LbmMode_>
                class SolverLBMGrid
                {
                };

        template<typename Tag_, typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
            class SolverLBMGrid<Tag_, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> : public SolverLBMGridBase
            {
                private:
                    unsigned long _parts;
                    PackedGridInfo<D2Q9> * _info;
                    PackedGridData<D2Q9, ResPrec_> * _data;
                    std::vector<PackedGridInfo<D2Q9> > _info_list;
                    std::vector<PackedGridData<D2Q9, ResPrec_> > _data_list;
                    std::vector<PackedGridFringe<D2Q9> > _fringe_list;
                    std::vector<honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> *> _solver_list;
                    std::vector<Ticket<tags::CPU::MultiCore> *> _tickets;

                public:
                    SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time):
                        _info(info),
                        _data(data)
                {
                    _parts = Configuration::instance()->get_value("mc::SolverLabsweGrid::patch_count", 4ul);
                    CONTEXT("When creating LABSWE solver:");
                    GridPartitioner<D2Q9, ResPrec_>::decompose(_parts, *_info, *_data, _info_list, _data_list, _fringe_list);

                    for(unsigned long i(0) ; i < _parts ; ++i)
                    {
                        _solver_list.push_back(new honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>(&_info_list[i], &_data_list[i], dx, dy, dt, rel_time));
                    }
                }

                    virtual ~SolverLBMGrid()
                    {
                        CONTEXT("When destroying LABSWE solver.");
                    }

                    void do_preprocessing()
                    {
                        CONTEXT("When performing LABSWE preprocessing.");

                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            _tickets.push_back(mc::ThreadPool::instance()->enqueue(
                                        bind(
                                            honei::GridPartitioner<D2Q9, ResPrec_>::recompose, &_info_list.at(i), &_data_list.at(i)),
                                            DispatchPolicy::on_core(i)));
                            tickets.push_back(_tickets.at(i));
                        }
                        tickets.wait();
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(mc::ThreadPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::do_preprocessing),
                                            *(_solver_list.at(i))
                                            ), DispatchPolicy::same_core_as(_tickets.at(i))));
                        }
                        tickets.wait();
                        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
                    }

                    void do_postprocessing()
                    {
                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(mc::ThreadPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::do_postprocessing),
                                            *(_solver_list.at(i))
                                            ), DispatchPolicy::same_core_as(_tickets.at(i))));
                        }
                        tickets.wait();
                        GridPartitioner<D2Q9, ResPrec_>::compose(*_info, *_data, _info_list, _data_list);
                    }

                    void solve()
                    {
                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(mc::ThreadPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::solve),
                                            *(_solver_list.at(i))
                                            ), DispatchPolicy::same_core_as(_tickets.at(i))));
                        }
                        tickets.wait();
                        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
                    }
            };
    }

    template<typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMGrid<tags::CPU::MultiCore, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> :
        public mc::SolverLBMGrid<tags::CPU::MultiCore, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>
        {
            public:
                SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time):
                    mc::SolverLBMGrid<tags::CPU::MultiCore, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>(info, data, dx, dy, dt, rel_time)
            {
            }
        };

    template<typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMGrid<tags::CPU::MultiCore::SSE, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> :
        public mc::SolverLBMGrid<tags::CPU::MultiCore::SSE, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>
        {
            public:
                SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time):
                    mc::SolverLBMGrid<tags::CPU::MultiCore::SSE, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>(info, data, dx, dy, dt, rel_time)
            {
            }
        };


    namespace
    {
        class UploadTask
        {
            private:
                void * device;
                void * address;
                unsigned long bytes;
            public:
                UploadTask(void * a, void * d, unsigned long b) :
                    device(d),
                    address(a),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_upload(address, device, bytes);
                }
        };

        class DownloadTask
        {
            private:
                void * device;
                void * address;
                unsigned long bytes;
            public:
                DownloadTask(void * d, void * a, unsigned long b) :
                    device(d),
                    address(a),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_download(device, address, bytes);
                }
        };

        class DetectTask
        {
            private:
                void ** device;
                void * address;
            public:
                DetectTask(void * a, void ** d) :
                    device(d),
                    address(a)
                {
                }

                void operator() ()
                {
                    *device = MemoryArbiter::instance()->lock(lm_read_only, tags::GPU::CUDA::memory_value, address, address, 0);
                    MemoryArbiter::instance()->unlock(lm_read_only, address);
                }
        };
    }
    namespace cuda
    {
        template<typename Tag_,
            typename Application_,
            typename ResPrec_,
            typename Force__,
            typename SourceScheme_,
            typename GridType_,
            typename LatticeType_,
            typename BoundaryType_,
            typename LbmMode_>
                class SolverLBMGrid
                {
                };

        template<typename Tag_, typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
            class SolverLBMGrid<Tag_, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> : public SolverLBMGridBase
            {
                private:
                    unsigned long _parts;
                    PackedGridInfo<D2Q9> * _info;
                    PackedGridData<D2Q9, ResPrec_> * _data;
                    std::vector<PackedGridInfo<D2Q9> > _info_list;
                    std::vector<PackedGridData<D2Q9, ResPrec_> > _data_list;
                    std::vector<PackedGridFringe<D2Q9> > _fringe_list;
                    std::vector<honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> *> _solver_list;

                public:
                    SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time):
                        _info(info),
                        _data(data)
                {
                    _parts = Configuration::instance()->get_value("cuda::SolverLabsweGrid::patch_count", 2ul);
                    CONTEXT("When creating LABSWE solver:");
                    GridPartitioner<D2Q9, ResPrec_>::decompose(_parts, *_info, *_data, _info_list, _data_list, _fringe_list);

                    for(unsigned long i(0) ; i < _parts ; ++i)
                    {
                        _solver_list.push_back(new honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>(&_info_list[i], &_data_list[i], dx, dy, dt, rel_time));
                    }
                }

                    virtual ~SolverLBMGrid()
                    {
                        CONTEXT("When destroying LABSWE solver.");
                    }

                    void do_preprocessing()
                    {
                        CONTEXT("When performing LABSWE preprocessing.");

                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::do_preprocessing),
                                            *(_solver_list.at(i))
                                            ), i));
                        }
                        tickets.wait();
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            _data_list[i].f_temp_1->lock(lm_read_and_write);
                            _data_list[i].f_temp_2->lock(lm_read_and_write);
                            _data_list[i].f_temp_3->lock(lm_read_and_write);
                            _data_list[i].f_temp_4->lock(lm_read_and_write);
                            _data_list[i].f_temp_5->lock(lm_read_and_write);
                            _data_list[i].f_temp_6->lock(lm_read_and_write);
                            _data_list[i].f_temp_7->lock(lm_read_and_write);
                            _data_list[i].f_temp_8->lock(lm_read_and_write);
                            _data_list[i].h->lock(lm_read_and_write);
                        }
                        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            _data_list[i].f_temp_1->unlock(lm_read_and_write);
                            _data_list[i].f_temp_2->unlock(lm_read_and_write);
                            _data_list[i].f_temp_3->unlock(lm_read_and_write);
                            _data_list[i].f_temp_4->unlock(lm_read_and_write);
                            _data_list[i].f_temp_5->unlock(lm_read_and_write);
                            _data_list[i].f_temp_6->unlock(lm_read_and_write);
                            _data_list[i].f_temp_7->unlock(lm_read_and_write);
                            _data_list[i].f_temp_8->unlock(lm_read_and_write);
                            _data_list[i].h->unlock(lm_read_and_write);
                        }
                    }

                    void do_postprocessing()
                    {
                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::do_postprocessing),
                                            *(_solver_list.at(i))
                                            ), i));
                        }
                        tickets.wait();
                        GridPartitioner<D2Q9, ResPrec_>::compose(*_info, *_data, _info_list, _data_list);
                    }

                    void solve()
                    {
                        TicketVector tickets;
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(
                                        bind(
                                            mem_fn(&honei::SolverLBMGrid<typename Tag_::DelegateTo, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>::solve),
                                            *(_solver_list.at(i))
                                            ), i));
                        }
                        tickets.wait();
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            //todo target nur einmal abrufen und dann speichern
                            void * target;
                            void * target1;
                            void * target2;
                            void * target3;
                            void * target4;
                            void * target5;
                            void * target6;
                            void * target7;
                            void * target8;
                            DetectTask detask(_data_list[i].h->address(), &target);
                            cuda::GPUPool::instance()->enqueue(detask, i)->wait();
                            DetectTask detask1(_data_list[i].f_temp_1->address(), &target1);
                            cuda::GPUPool::instance()->enqueue(detask1, i)->wait();
                            DetectTask detask2(_data_list[i].f_temp_2->address(), &target2);
                            cuda::GPUPool::instance()->enqueue(detask2, i)->wait();
                            DetectTask detask3(_data_list[i].f_temp_3->address(), &target3);
                            cuda::GPUPool::instance()->enqueue(detask3, i)->wait();
                            DetectTask detask4(_data_list[i].f_temp_4->address(), &target4);
                            cuda::GPUPool::instance()->enqueue(detask4, i)->wait();
                            DetectTask detask5(_data_list[i].f_temp_5->address(), &target5);
                            cuda::GPUPool::instance()->enqueue(detask5, i)->wait();
                            DetectTask detask6(_data_list[i].f_temp_6->address(), &target6);
                            cuda::GPUPool::instance()->enqueue(detask6, i)->wait();
                            DetectTask detask7(_data_list[i].f_temp_7->address(), &target7);
                            cuda::GPUPool::instance()->enqueue(detask7, i)->wait();
                            DetectTask detask8(_data_list[i].f_temp_8->address(), &target8);
                            cuda::GPUPool::instance()->enqueue(detask8, i)->wait();

                            unsigned long offset(_info_list[i].offset);
                            unsigned long f1_offset((*_fringe_list[i].dir_index_1)[0]);
                            unsigned long f1_size((*_fringe_list[i].dir_index_1)[_fringe_list[i].dir_index_1->size()-1] - f1_offset);
                            DownloadTask dtask1((void *)((ResPrec_*)target1 + f1_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_1->address() + f1_offset - offset), f1_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask1, i));
                            unsigned long f2_offset((*_fringe_list[i].dir_index_2)[0]);
                            unsigned long f2_size((*_fringe_list[i].dir_index_2)[_fringe_list[i].dir_index_2->size()-1] - f2_offset);
                            DownloadTask dtask2((void *)((ResPrec_*)target2 + f2_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_2->address() + f2_offset - offset), f2_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask2, i));
                            unsigned long f3_offset((*_fringe_list[i].dir_index_3)[0]);
                            unsigned long f3_size((*_fringe_list[i].dir_index_3)[_fringe_list[i].dir_index_3->size()-1] - f3_offset);
                            DownloadTask dtask3((void *)((ResPrec_*)target3 + f3_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_3->address() + f3_offset - offset), f3_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask3, i));
                            unsigned long f4_offset((*_fringe_list[i].dir_index_4)[0]);
                            unsigned long f4_size((*_fringe_list[i].dir_index_4)[_fringe_list[i].dir_index_4->size()-1] - f4_offset);
                            DownloadTask dtask4((void *)((ResPrec_*)target4 + f4_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_4->address() + f4_offset - offset), f4_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask4, i));
                            unsigned long f5_offset((*_fringe_list[i].dir_index_5)[0]);
                            unsigned long f5_size((*_fringe_list[i].dir_index_5)[_fringe_list[i].dir_index_5->size()-1] - f5_offset);
                            DownloadTask dtask5((void *)((ResPrec_*)target5 + f5_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_5->address() + f5_offset - offset), f5_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask5, i));
                            unsigned long f6_offset((*_fringe_list[i].dir_index_6)[0]);
                            unsigned long f6_size((*_fringe_list[i].dir_index_6)[_fringe_list[i].dir_index_6->size()-1] - f6_offset);
                            DownloadTask dtask6((void *)((ResPrec_*)target6 + f6_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_6->address() + f6_offset - offset), f6_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask6, i));
                            unsigned long f7_offset((*_fringe_list[i].dir_index_7)[0]);
                            unsigned long f7_size((*_fringe_list[i].dir_index_7)[_fringe_list[i].dir_index_7->size()-1] - f7_offset);
                            DownloadTask dtask7((void *)((ResPrec_*)target7 + f7_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_7->address() + f7_offset - offset), f7_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask7, i));
                            unsigned long f8_offset((*_fringe_list[i].dir_index_8)[0]);
                            unsigned long f8_size((*_fringe_list[i].dir_index_8)[_fringe_list[i].dir_index_8->size()-1] - f8_offset);
                            DownloadTask dtask8((void *)((ResPrec_*)target8 + f8_offset - offset), (void *)((ResPrec_*)_data_list[i].f_temp_8->address() + f8_offset - offset), f8_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask8, i));

                            for (unsigned long j(0) ; j < _fringe_list[i].external_h_index->size() / 2 ; ++j)
                            {
                                unsigned long h_offset((*_fringe_list[i].external_h_index)[j * 2]);
                                unsigned long h_size((*_fringe_list[i].external_h_index)[j * 2 + 1] - h_offset);
                                DownloadTask dtask((void *)((ResPrec_*)target + h_offset - offset), (void *)((ResPrec_*)_data_list[i].h->address() + h_offset - offset), h_size * sizeof(ResPrec_));
                                tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask, i));
                            }

                        }
                        tickets.wait();
                        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
                        for (unsigned long i(0) ; i < _parts ; ++i)
                        {
                            void * target;
                            void * target1;
                            void * target2;
                            void * target3;
                            void * target4;
                            void * target5;
                            void * target6;
                            void * target7;
                            void * target8;
                            DetectTask detask(_data_list[i].h->address(), &target);
                            cuda::GPUPool::instance()->enqueue(detask, i)->wait();
                            DetectTask detask1(_data_list[i].f_temp_1->address(), &target1);
                            cuda::GPUPool::instance()->enqueue(detask1, i)->wait();
                            DetectTask detask2(_data_list[i].f_temp_2->address(), &target2);
                            cuda::GPUPool::instance()->enqueue(detask2, i)->wait();
                            DetectTask detask3(_data_list[i].f_temp_3->address(), &target3);
                            cuda::GPUPool::instance()->enqueue(detask3, i)->wait();
                            DetectTask detask4(_data_list[i].f_temp_4->address(), &target4);
                            cuda::GPUPool::instance()->enqueue(detask4, i)->wait();
                            DetectTask detask5(_data_list[i].f_temp_5->address(), &target5);
                            cuda::GPUPool::instance()->enqueue(detask5, i)->wait();
                            DetectTask detask6(_data_list[i].f_temp_6->address(), &target6);
                            cuda::GPUPool::instance()->enqueue(detask6, i)->wait();
                            DetectTask detask7(_data_list[i].f_temp_7->address(), &target7);
                            cuda::GPUPool::instance()->enqueue(detask7, i)->wait();
                            DetectTask detask8(_data_list[i].f_temp_8->address(), &target8);
                            cuda::GPUPool::instance()->enqueue(detask8, i)->wait();

                            unsigned long offset(_info_list[i].offset);
                            unsigned long f1_offset((*_fringe_list[i].external_dir_index_1)[0]);
                            unsigned long f1_size((*_fringe_list[i].external_dir_index_1)[_fringe_list[i].external_dir_index_1->size()-1] - f1_offset);
                            UploadTask utask1((void *)((ResPrec_*)_data_list[i].f_temp_1->address() + f1_offset - offset), (void *)((ResPrec_*)target1 + f1_offset - offset), f1_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask1, i));
                            unsigned long f2_offset((*_fringe_list[i].external_dir_index_2)[0]);
                            unsigned long f2_size((*_fringe_list[i].external_dir_index_2)[_fringe_list[i].external_dir_index_2->size()-1] - f2_offset);
                            UploadTask utask2((void *)((ResPrec_*)_data_list[i].f_temp_2->address() + f2_offset - offset), (void *)((ResPrec_*)target2 + f2_offset - offset), f2_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask2, i));
                            unsigned long f3_offset((*_fringe_list[i].external_dir_index_3)[0]);
                            unsigned long f3_size((*_fringe_list[i].external_dir_index_3)[_fringe_list[i].external_dir_index_3->size()-1] - f3_offset);
                            UploadTask utask3((void *)((ResPrec_*)_data_list[i].f_temp_3->address() + f3_offset - offset), (void *)((ResPrec_*)target3 + f3_offset - offset), f3_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask3, i));
                            unsigned long f4_offset((*_fringe_list[i].external_dir_index_4)[0]);
                            unsigned long f4_size((*_fringe_list[i].external_dir_index_4)[_fringe_list[i].external_dir_index_4->size()-1] - f4_offset);
                            UploadTask utask4((void *)((ResPrec_*)_data_list[i].f_temp_4->address() + f4_offset - offset), (void *)((ResPrec_*)target4 + f4_offset - offset), f4_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask4, i));
                            unsigned long f5_offset((*_fringe_list[i].external_dir_index_5)[0]);
                            unsigned long f5_size((*_fringe_list[i].external_dir_index_5)[_fringe_list[i].external_dir_index_5->size()-1] - f5_offset);
                            UploadTask utask5((void *)((ResPrec_*)_data_list[i].f_temp_5->address() + f5_offset - offset), (void *)((ResPrec_*)target5 + f5_offset - offset), f5_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask5, i));
                            unsigned long f6_offset((*_fringe_list[i].external_dir_index_6)[0]);
                            unsigned long f6_size((*_fringe_list[i].external_dir_index_6)[_fringe_list[i].external_dir_index_6->size()-1] - f6_offset);
                            UploadTask utask6((void *)((ResPrec_*)_data_list[i].f_temp_6->address() + f6_offset - offset), (void *)((ResPrec_*)target6 + f6_offset - offset), f6_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask6, i));
                            unsigned long f7_offset((*_fringe_list[i].external_dir_index_7)[0]);
                            unsigned long f7_size((*_fringe_list[i].external_dir_index_7)[_fringe_list[i].external_dir_index_7->size()-1] - f7_offset);
                            UploadTask utask7((void *)((ResPrec_*)_data_list[i].f_temp_7->address() + f7_offset - offset), (void *)((ResPrec_*)target7 + f7_offset - offset), f7_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask7, i));
                            unsigned long f8_offset((*_fringe_list[i].external_dir_index_8)[0]);
                            unsigned long f8_size((*_fringe_list[i].external_dir_index_8)[_fringe_list[i].external_dir_index_8->size()-1] - f8_offset);
                            UploadTask utask8((void *)((ResPrec_*)_data_list[i].f_temp_8->address() + f8_offset - offset), (void *)((ResPrec_*)target8 + f8_offset - offset), f8_size * sizeof(ResPrec_));
                            tickets.push_back(cuda::GPUPool::instance()->enqueue(utask8, i));

                            for (unsigned long j(0) ; j < _fringe_list[i].h_index->size() / 2 ; ++j)
                            {
                                unsigned long h_offset((*_fringe_list[i].h_index)[j * 2]);
                                unsigned long h_size((*_fringe_list[i].h_index)[j * 2 + 1] - h_offset);
                                UploadTask utask((void *)((ResPrec_*)_data_list[i].h->address() + h_offset - offset), (void *)((ResPrec_*)target + h_offset - offset), h_size * sizeof(ResPrec_));
                                tickets.push_back(cuda::GPUPool::instance()->enqueue(utask, i));
                            }
                        }
                        tickets.wait();
                    }
            };
    }

    template<typename Application_, typename ResPrec_, typename Force_, typename SourceScheme_, typename LbmMode_>
        class SolverLBMGrid<tags::GPU::MultiCore::CUDA, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_> :
        public cuda::SolverLBMGrid<tags::GPU::MultiCore::CUDA, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>
        {
            public:
                SolverLBMGrid(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, ResPrec_> * data, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, ResPrec_ rel_time):
                    cuda::SolverLBMGrid<tags::GPU::MultiCore::CUDA, Application_, ResPrec_, Force_, SourceScheme_, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, LbmMode_>(info, data, dx, dy, dt, rel_time)
            {
            }
        };
}
#endif
