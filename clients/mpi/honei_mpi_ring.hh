/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
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

#ifndef MPI_GUARD_HONEI_MPI_HH
#define MPU_GUARD_HONEI_MPI_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/swe/volume.hh>
#include <unittest/unittest.hh>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/grid_partitioner.hh>
#include <honei/util/time_stamp.hh>
#include <honei/lbm/scenario_collection.hh>

#include <iostream>
#include <vector>



#include <mpi.h>

namespace honei
{
    template <typename Tag_, typename DataType_>
    class MPIRingSolver
    {
        private:
            int _numprocs;
            int _myid;

        public:
            MPIRingSolver(int argc, char **argv)
            {
                unsigned long gridsize(atoi(argv[1]));
                mpi::mpi_init(&argc, &argv);
                mpi::mpi_comm_size(&_numprocs);
                mpi::mpi_comm_rank(&_myid);

                if (_myid == 0)
                {
                    _master(gridsize);
                }
                else
                {
                    _slave();
                }
            }

            ~MPIRingSolver()
            {
                mpi::mpi_finalize();
            }

        private:
            void _master(unsigned long gridsize)
            {
                std::cout<<"Ring LBM Solver with " << _numprocs << " nodes:" << std::endl;
                unsigned long timesteps(100);
                Grid<D2Q9, DataType_> grid;
                ScenarioCollection::get_scenario(0, gridsize, gridsize, grid);
                std::cout << "Solving: " << grid.description << std::endl;

                PackedGridData<D2Q9, DataType_>  data;
                PackedGridInfo<D2Q9> info;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
                std::vector<PackedGridInfo<D2Q9> > info_list;
                std::vector<PackedGridData<D2Q9, DataType_> > data_list;
                std::vector<PackedGridFringe<D2Q9> > fringe_list;
                GridPartitioner<D2Q9, DataType_>::decompose(_numprocs - 1, info, data, info_list, data_list, fringe_list);

                mpi::mpi_bcast(&timesteps, 1, 0);
                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _send_info(target, info_list[target - 1]);
                    _send_data(target, data_list[target - 1]);
                    _send_fringe(target, fringe_list[target - 1]);
                }


                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _recv_slave_sync(target, info_list[target - 1], data_list[target - 1], fringe_list[target - 1]);
                }

                GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list);

                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _send_slave_sync(target, info_list[target - 1], data_list[target - 1], fringe_list[target - 1]);
                }

                //GridPartitioner<D2Q9, DataType_>::compose(info, data, info_list, data_list);
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                //PostProcessing<output_types::GNUPLOT>::value(*grid.h, 1, grid.h->columns(), grid.h->rows(), 101);
                // Preproc finished
                // start timesteps
                TimeStamp at, bt;
                at.take();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    //here are the solvers solving...
                    //and finished
                    /*for (signed long target(1) ; target < _numprocs ; ++target)
                    {
                        _recv_slave_sync(target, info_list[target - 1], data_list[target - 1], fringe_list[target - 1]);
                    }

                    GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list);

                    for (signed long target(1) ; target < _numprocs ; ++target)
                    {
                        _send_slave_sync(target, info_list[target - 1], data_list[target - 1], fringe_list[target - 1]);
                    }*/
                }
                MPI_Barrier(MPI_COMM_WORLD);
                bt.take();
                std::cout<<"Gridsize: "<<grid.h->rows()<<" x "<<grid.h->columns()<<std::endl;
                std::cout<<"Timesteps: " << timesteps << " TOE: "<<bt.total() - at.total()<<std::endl;
                std::cout<<"MLUPS: "<< (double(grid.h->rows()) * double(grid.h->columns()) * double(timesteps)) / (1e6 * (bt.total() - at.total())) <<std::endl;

                // generate output
                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _recv_full_sync(target, data_list[target - 1]);
                }
                GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list);
                GridPartitioner<D2Q9, DataType_>::compose(info, data, info_list, data_list);
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                //std::cout<<*grid.h;
                //PostProcessing<output_types::GNUPLOT>::value(*grid.h, 1, grid.h->columns(), grid.h->rows(),

                Grid<D2Q9, DataType_> grid_ref;
                ScenarioCollection::get_scenario(0, gridsize, gridsize, grid_ref);

                PackedGridData<D2Q9, DataType_>  data_ref;
                PackedGridInfo<D2Q9> info_ref;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_ref, info_ref, data_ref);
                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info_ref, &data_ref, 0.01, 0.01, 0.01, 1.1);
                solver.do_preprocessing();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver.solve();
                }
                for (unsigned long i(0) ; i < data.h->size() ; ++i)
                {
                    if (fabs((*data.h)[i] - (*data_ref.h)[i]) > 0.0001)
                        std::cout<<(*data.h)[i]<<" "<<(*data_ref.h)[i]<<std::endl;
                }
            }

            void _slave()
            {
                PackedGridData<D2Q9, DataType_> data;
                PackedGridInfo<D2Q9> info;
                PackedGridFringe<D2Q9> fringe;
                unsigned long timesteps;

                mpi::mpi_bcast(&timesteps, 1, 0);

                _recv_info(info);
                _recv_data(data);
                _recv_fringe(fringe);

                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, 0.01, 0.01, 0.01, 1.1);

                solver.do_preprocessing();

                _send_master_sync(0, info, data, fringe);
                _recv_master_sync(0, info, data, fringe);

                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver.solve();

                    //_send_master_sync(0, info, data, fringe);
                    //_recv_master_sync(0, info, data, fringe);

                    _circle_sync(info, data, fringe);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                _send_full_sync(0, data);
            }

            void _send_info(unsigned long target, PackedGridInfo<D2Q9> & info)
            {
                unsigned long limits_size(info.limits->size());
                unsigned long dir_index_1_size(info.dir_index_1->size());
                unsigned long dir_1_size(info.dir_1->size());
                unsigned long dir_index_2_size(info.dir_index_2->size());
                unsigned long dir_2_size(info.dir_2->size());
                unsigned long dir_index_3_size(info.dir_index_3->size());
                unsigned long dir_3_size(info.dir_3->size());
                unsigned long dir_index_4_size(info.dir_index_4->size());
                unsigned long dir_4_size(info.dir_4->size());
                unsigned long dir_index_5_size(info.dir_index_5->size());
                unsigned long dir_5_size(info.dir_5->size());
                unsigned long dir_index_6_size(info.dir_index_6->size());
                unsigned long dir_6_size(info.dir_6->size());
                unsigned long dir_index_7_size(info.dir_index_7->size());
                unsigned long dir_7_size(info.dir_7->size());
                unsigned long dir_index_8_size(info.dir_index_8->size());
                unsigned long dir_8_size(info.dir_8->size());

                mpi::mpi_send(&info.offset, 1, target, _myid);
                mpi::mpi_send(&limits_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_1_size, 1, target, _myid);
                mpi::mpi_send(&dir_1_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_2_size, 1, target, _myid);
                mpi::mpi_send(&dir_2_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_3_size, 1, target, _myid);
                mpi::mpi_send(&dir_3_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_4_size, 1, target, _myid);
                mpi::mpi_send(&dir_4_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_5_size, 1, target, _myid);
                mpi::mpi_send(&dir_5_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_6_size, 1, target, _myid);
                mpi::mpi_send(&dir_6_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_7_size, 1, target, _myid);
                mpi::mpi_send(&dir_7_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_8_size, 1, target, _myid);
                mpi::mpi_send(&dir_8_size, 1, target, _myid);

                mpi::mpi_send(info.limits->elements(), info.limits->size(), target, _myid);
                mpi::mpi_send(info.types->elements(), info.types->size(), target, _myid);
                mpi::mpi_send(info.dir_index_1->elements(), info.dir_index_1->size(), target, _myid);
                mpi::mpi_send(info.dir_1->elements(), info.dir_1->size(), target, _myid);
                mpi::mpi_send(info.dir_index_2->elements(), info.dir_index_2->size(), target, _myid);
                mpi::mpi_send(info.dir_2->elements(), info.dir_2->size(), target, _myid);
                mpi::mpi_send(info.dir_index_3->elements(), info.dir_index_3->size(), target, _myid);
                mpi::mpi_send(info.dir_3->elements(), info.dir_3->size(), target, _myid);
                mpi::mpi_send(info.dir_index_4->elements(), info.dir_index_4->size(), target, _myid);
                mpi::mpi_send(info.dir_4->elements(), info.dir_4->size(), target, _myid);
                mpi::mpi_send(info.dir_index_5->elements(), info.dir_index_5->size(), target, _myid);
                mpi::mpi_send(info.dir_5->elements(), info.dir_5->size(), target, _myid);
                mpi::mpi_send(info.dir_index_6->elements(), info.dir_index_6->size(), target, _myid);
                mpi::mpi_send(info.dir_6->elements(), info.dir_6->size(), target, _myid);
                mpi::mpi_send(info.dir_index_7->elements(), info.dir_index_7->size(), target, _myid);
                mpi::mpi_send(info.dir_7->elements(), info.dir_7->size(), target, _myid);
                mpi::mpi_send(info.dir_index_8->elements(), info.dir_index_8->size(), target, _myid);
                mpi::mpi_send(info.dir_8->elements(), info.dir_8->size(), target, _myid);
            }

            void _recv_info(PackedGridInfo<D2Q9> & info)
            {
                unsigned long limits_size;
                unsigned long dir_index_1_size;
                unsigned long dir_1_size;
                unsigned long dir_index_2_size;
                unsigned long dir_2_size;
                unsigned long dir_index_3_size;
                unsigned long dir_3_size;
                unsigned long dir_index_4_size;
                unsigned long dir_4_size;
                unsigned long dir_index_5_size;
                unsigned long dir_5_size;
                unsigned long dir_index_6_size;
                unsigned long dir_6_size;
                unsigned long dir_index_7_size;
                unsigned long dir_7_size;
                unsigned long dir_index_8_size;
                unsigned long dir_8_size;

                mpi::mpi_recv(&info.offset, 1, 0, 0);
                mpi::mpi_recv(&limits_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_1_size, 1, 0, 0);
                mpi::mpi_recv(&dir_1_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_2_size, 1, 0, 0);
                mpi::mpi_recv(&dir_2_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_3_size, 1, 0, 0);
                mpi::mpi_recv(&dir_3_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_4_size, 1, 0, 0);
                mpi::mpi_recv(&dir_4_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_5_size, 1, 0, 0);
                mpi::mpi_recv(&dir_5_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_6_size, 1, 0, 0);
                mpi::mpi_recv(&dir_6_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_7_size, 1, 0, 0);
                mpi::mpi_recv(&dir_7_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_8_size, 1, 0, 0);
                mpi::mpi_recv(&dir_8_size, 1, 0, 0);

                info.limits = new DenseVector<unsigned long>(limits_size);
                info.types = new DenseVector<unsigned long>(limits_size);
                info.dir_index_1 = new DenseVector<unsigned long>(dir_index_1_size);
                info.dir_1 = new DenseVector<unsigned long>(dir_1_size);
                info.dir_index_2 = new DenseVector<unsigned long>(dir_index_2_size);
                info.dir_2 = new DenseVector<unsigned long>(dir_2_size);
                info.dir_index_3 = new DenseVector<unsigned long>(dir_index_3_size);
                info.dir_3 = new DenseVector<unsigned long>(dir_3_size);
                info.dir_index_4 = new DenseVector<unsigned long>(dir_index_4_size);
                info.dir_4 = new DenseVector<unsigned long>(dir_4_size);
                info.dir_index_5 = new DenseVector<unsigned long>(dir_index_5_size);
                info.dir_5 = new DenseVector<unsigned long>(dir_5_size);
                info.dir_index_6 = new DenseVector<unsigned long>(dir_index_6_size);
                info.dir_6 = new DenseVector<unsigned long>(dir_6_size);
                info.dir_index_7 = new DenseVector<unsigned long>(dir_index_7_size);
                info.dir_7 = new DenseVector<unsigned long>(dir_7_size);
                info.dir_index_8 = new DenseVector<unsigned long>(dir_index_8_size);
                info.dir_8 = new DenseVector<unsigned long>(dir_8_size);

                mpi::mpi_recv(info.limits->elements(), info.limits->size(), 0, 0);
                mpi::mpi_recv(info.types->elements(), info.types->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_1->elements(), info.dir_index_1->size(), 0, 0);
                mpi::mpi_recv(info.dir_1->elements(), info.dir_1->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_2->elements(), info.dir_index_2->size(), 0, 0);
                mpi::mpi_recv(info.dir_2->elements(), info.dir_2->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_3->elements(), info.dir_index_3->size(), 0, 0);
                mpi::mpi_recv(info.dir_3->elements(), info.dir_3->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_4->elements(), info.dir_index_4->size(), 0, 0);
                mpi::mpi_recv(info.dir_4->elements(), info.dir_4->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_5->elements(), info.dir_index_5->size(), 0, 0);
                mpi::mpi_recv(info.dir_5->elements(), info.dir_5->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_6->elements(), info.dir_index_6->size(), 0, 0);
                mpi::mpi_recv(info.dir_6->elements(), info.dir_6->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_7->elements(), info.dir_index_7->size(), 0, 0);
                mpi::mpi_recv(info.dir_7->elements(), info.dir_7->size(), 0, 0);
                mpi::mpi_recv(info.dir_index_8->elements(), info.dir_index_8->size(), 0, 0);
                mpi::mpi_recv(info.dir_8->elements(), info.dir_8->size(), 0, 0);
            }

            void _send_data(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                unsigned long h_size(data.h->size());
                unsigned long dist_size(data.distribution_x->size());

                mpi::mpi_send(&h_size, 1, target, _myid);
                mpi::mpi_send(&dist_size, 1, target, _myid);

                mpi::mpi_send(data.h->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.u->elements(), data.u->size(), target, _myid);
                mpi::mpi_send(data.v->elements(), data.v->size(), target, _myid);
                mpi::mpi_send(data.b->elements(), data.b->size(), target, _myid);
            }

            void _recv_data(PackedGridData<D2Q9, DataType_> & data)
            {
                unsigned long h_size;
                unsigned long dist_size;

                mpi::mpi_recv(&h_size, 1, 0, 0);
                mpi::mpi_recv(&dist_size, 1, 0, 0);

                data.h = new DenseVector<DataType_>(h_size);
                data.u = new DenseVector<DataType_>(h_size);
                data.v = new DenseVector<DataType_>(h_size);
                data.b = new DenseVector<DataType_>(h_size);
                data.temp = new DenseVector<DataType_>(h_size);

                mpi::mpi_recv(data.h->elements(), data.h->size(), 0, 0);
                mpi::mpi_recv(data.u->elements(), data.u->size(), 0, 0);
                mpi::mpi_recv(data.v->elements(), data.v->size(), 0, 0);
                mpi::mpi_recv(data.b->elements(), data.b->size(), 0, 0);

                data.f_0 = new DenseVector<DataType_>(h_size);
                data.f_1 = new DenseVector<DataType_>(h_size);
                data.f_2 = new DenseVector<DataType_>(h_size);
                data.f_3 = new DenseVector<DataType_>(h_size);
                data.f_4 = new DenseVector<DataType_>(h_size);
                data.f_5 = new DenseVector<DataType_>(h_size);
                data.f_6 = new DenseVector<DataType_>(h_size);
                data.f_7 = new DenseVector<DataType_>(h_size);
                data.f_8 = new DenseVector<DataType_>(h_size);

                data.f_eq_0 = new DenseVector<DataType_>(h_size);
                data.f_eq_1 = new DenseVector<DataType_>(h_size);
                data.f_eq_2 = new DenseVector<DataType_>(h_size);
                data.f_eq_3 = new DenseVector<DataType_>(h_size);
                data.f_eq_4 = new DenseVector<DataType_>(h_size);
                data.f_eq_5 = new DenseVector<DataType_>(h_size);
                data.f_eq_6 = new DenseVector<DataType_>(h_size);
                data.f_eq_7 = new DenseVector<DataType_>(h_size);
                data.f_eq_8 = new DenseVector<DataType_>(h_size);

                data.f_temp_0 = new DenseVector<DataType_>(h_size);
                data.f_temp_1 = new DenseVector<DataType_>(h_size);
                data.f_temp_2 = new DenseVector<DataType_>(h_size);
                data.f_temp_3 = new DenseVector<DataType_>(h_size);
                data.f_temp_4 = new DenseVector<DataType_>(h_size);
                data.f_temp_5 = new DenseVector<DataType_>(h_size);
                data.f_temp_6 = new DenseVector<DataType_>(h_size);
                data.f_temp_7 = new DenseVector<DataType_>(h_size);
                data.f_temp_8 = new DenseVector<DataType_>(h_size);

                data.distribution_x = new DenseVector<DataType_>(dist_size);
                data.distribution_y = new DenseVector<DataType_>(dist_size);
            }

            void _send_fringe(unsigned long target, PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long h_targets_size(fringe.h_targets->size());
                unsigned long h_index_size(fringe.h_index->size());
                unsigned long external_h_index_size(fringe.external_h_index->size());
                unsigned long dir_targets_1_size(fringe.dir_targets_1->size());
                unsigned long dir_index_1_size(fringe.dir_index_1->size());
                unsigned long dir_targets_2_size(fringe.dir_targets_2->size());
                unsigned long dir_index_2_size(fringe.dir_index_2->size());
                unsigned long dir_targets_3_size(fringe.dir_targets_3->size());
                unsigned long dir_index_3_size(fringe.dir_index_3->size());
                unsigned long dir_targets_4_size(fringe.dir_targets_4->size());
                unsigned long dir_index_4_size(fringe.dir_index_4->size());
                unsigned long dir_targets_5_size(fringe.dir_targets_5->size());
                unsigned long dir_index_5_size(fringe.dir_index_5->size());
                unsigned long dir_targets_6_size(fringe.dir_targets_6->size());
                unsigned long dir_index_6_size(fringe.dir_index_6->size());
                unsigned long dir_targets_7_size(fringe.dir_targets_7->size());
                unsigned long dir_index_7_size(fringe.dir_index_7->size());
                unsigned long dir_targets_8_size(fringe.dir_targets_8->size());
                unsigned long dir_index_8_size(fringe.dir_index_8->size());
                unsigned long external_dir_index_1_size(fringe.external_dir_index_1->size());
                unsigned long external_dir_index_2_size(fringe.external_dir_index_2->size());
                unsigned long external_dir_index_3_size(fringe.external_dir_index_3->size());
                unsigned long external_dir_index_4_size(fringe.external_dir_index_4->size());
                unsigned long external_dir_index_5_size(fringe.external_dir_index_5->size());
                unsigned long external_dir_index_6_size(fringe.external_dir_index_6->size());
                unsigned long external_dir_index_7_size(fringe.external_dir_index_7->size());
                unsigned long external_dir_index_8_size(fringe.external_dir_index_8->size());
                unsigned long external_h_targets_size(fringe.external_h_targets->size());
                unsigned long external_dir_targets_1_size(fringe.external_dir_targets_1->size());
                unsigned long external_dir_targets_2_size(fringe.external_dir_targets_2->size());
                unsigned long external_dir_targets_3_size(fringe.external_dir_targets_3->size());
                unsigned long external_dir_targets_4_size(fringe.external_dir_targets_4->size());
                unsigned long external_dir_targets_5_size(fringe.external_dir_targets_5->size());
                unsigned long external_dir_targets_6_size(fringe.external_dir_targets_6->size());
                unsigned long external_dir_targets_7_size(fringe.external_dir_targets_7->size());
                unsigned long external_dir_targets_8_size(fringe.external_dir_targets_8->size());

                mpi::mpi_send(&h_targets_size, 1, target, _myid);
                mpi::mpi_send(&h_index_size, 1, target, _myid);
                mpi::mpi_send(&external_h_index_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_1_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_1_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_2_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_2_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_3_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_3_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_4_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_4_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_5_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_5_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_6_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_6_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_7_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_7_size, 1, target, _myid);
                mpi::mpi_send(&dir_targets_8_size, 1, target, _myid);
                mpi::mpi_send(&dir_index_8_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_1_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_2_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_3_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_4_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_5_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_6_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_7_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_index_8_size, 1, target, _myid);
                mpi::mpi_send(&external_h_targets_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_1_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_2_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_3_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_4_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_5_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_6_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_7_size, 1, target, _myid);
                mpi::mpi_send(&external_dir_targets_8_size, 1, target, _myid);

                mpi::mpi_send(fringe.h_targets->elements(), fringe.h_targets->size(), target, _myid);
                mpi::mpi_send(fringe.h_index->elements(), fringe.h_index->size(), target, _myid);
                mpi::mpi_send(fringe.external_h_index->elements(), fringe.external_h_index->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_1->elements(), fringe.dir_targets_1->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_1->elements(), fringe.dir_index_1->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_2->elements(), fringe.dir_targets_2->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_2->elements(), fringe.dir_index_2->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_3->elements(), fringe.dir_targets_3->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_3->elements(), fringe.dir_index_3->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_4->elements(), fringe.dir_targets_4->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_4->elements(), fringe.dir_index_4->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_5->elements(), fringe.dir_targets_5->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_5->elements(), fringe.dir_index_5->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_6->elements(), fringe.dir_targets_6->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_6->elements(), fringe.dir_index_6->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_7->elements(), fringe.dir_targets_7->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_7->elements(), fringe.dir_index_7->size(), target, _myid);
                mpi::mpi_send(fringe.dir_targets_8->elements(), fringe.dir_targets_8->size(), target, _myid);
                mpi::mpi_send(fringe.dir_index_8->elements(), fringe.dir_index_8->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_1->elements(), fringe.external_dir_index_1->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_2->elements(), fringe.external_dir_index_2->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_3->elements(), fringe.external_dir_index_3->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_4->elements(), fringe.external_dir_index_4->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_5->elements(), fringe.external_dir_index_5->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_6->elements(), fringe.external_dir_index_6->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_7->elements(), fringe.external_dir_index_7->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_index_8->elements(), fringe.external_dir_index_8->size(), target, _myid);
                mpi::mpi_send(fringe.external_h_targets->elements(), fringe.external_h_targets->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_1->elements(), fringe.external_dir_targets_1->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_2->elements(), fringe.external_dir_targets_2->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_3->elements(), fringe.external_dir_targets_3->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_4->elements(), fringe.external_dir_targets_4->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_5->elements(), fringe.external_dir_targets_5->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_6->elements(), fringe.external_dir_targets_6->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_7->elements(), fringe.external_dir_targets_7->size(), target, _myid);
                mpi::mpi_send(fringe.external_dir_targets_8->elements(), fringe.external_dir_targets_8->size(), target, _myid);
            }

            void _recv_fringe(PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long h_targets_size;
                unsigned long h_index_size;
                unsigned long external_h_index_size;
                unsigned long dir_targets_1_size;
                unsigned long dir_index_1_size;
                unsigned long dir_targets_2_size;
                unsigned long dir_index_2_size;
                unsigned long dir_targets_3_size;
                unsigned long dir_index_3_size;
                unsigned long dir_targets_4_size;
                unsigned long dir_index_4_size;
                unsigned long dir_targets_5_size;
                unsigned long dir_index_5_size;
                unsigned long dir_targets_6_size;
                unsigned long dir_index_6_size;
                unsigned long dir_targets_7_size;
                unsigned long dir_index_7_size;
                unsigned long dir_targets_8_size;
                unsigned long dir_index_8_size;
                unsigned long external_dir_index_1_size;
                unsigned long external_dir_index_2_size;
                unsigned long external_dir_index_3_size;
                unsigned long external_dir_index_4_size;
                unsigned long external_dir_index_5_size;
                unsigned long external_dir_index_6_size;
                unsigned long external_dir_index_7_size;
                unsigned long external_dir_index_8_size;
                unsigned long external_h_targets_size;
                unsigned long external_dir_targets_1_size;
                unsigned long external_dir_targets_2_size;
                unsigned long external_dir_targets_3_size;
                unsigned long external_dir_targets_4_size;
                unsigned long external_dir_targets_5_size;
                unsigned long external_dir_targets_6_size;
                unsigned long external_dir_targets_7_size;
                unsigned long external_dir_targets_8_size;

                mpi::mpi_recv(&h_targets_size, 1, 0, 0);
                mpi::mpi_recv(&h_index_size, 1, 0, 0);
                mpi::mpi_recv(&external_h_index_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_1_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_1_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_2_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_2_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_3_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_3_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_4_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_4_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_5_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_5_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_6_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_6_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_7_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_7_size, 1, 0, 0);
                mpi::mpi_recv(&dir_targets_8_size, 1, 0, 0);
                mpi::mpi_recv(&dir_index_8_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_1_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_2_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_3_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_4_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_5_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_6_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_7_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_index_8_size, 1, 0, 0);
                mpi::mpi_recv(&external_h_targets_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_1_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_2_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_3_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_4_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_5_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_6_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_7_size, 1, 0, 0);
                mpi::mpi_recv(&external_dir_targets_8_size, 1, 0, 0);

                fringe.h_targets = new DenseVector<unsigned long>(h_targets_size);
                fringe.h_index = new DenseVector<unsigned long>(h_index_size);
                fringe.external_h_index = new DenseVector<unsigned long>(external_h_index_size);
                fringe.dir_targets_1 = new DenseVector<unsigned long>(dir_targets_1_size);
                fringe.dir_index_1 = new DenseVector<unsigned long>(dir_index_1_size);
                fringe.dir_targets_2 = new DenseVector<unsigned long>(dir_targets_2_size);
                fringe.dir_index_2 = new DenseVector<unsigned long>(dir_index_2_size);
                fringe.dir_targets_3 = new DenseVector<unsigned long>(dir_targets_3_size);
                fringe.dir_index_3 = new DenseVector<unsigned long>(dir_index_3_size);
                fringe.dir_targets_4 = new DenseVector<unsigned long>(dir_targets_4_size);
                fringe.dir_index_4 = new DenseVector<unsigned long>(dir_index_4_size);
                fringe.dir_targets_5 = new DenseVector<unsigned long>(dir_targets_5_size);
                fringe.dir_index_5 = new DenseVector<unsigned long>(dir_index_5_size);
                fringe.dir_targets_6 = new DenseVector<unsigned long>(dir_targets_6_size);
                fringe.dir_index_6 = new DenseVector<unsigned long>(dir_index_6_size);
                fringe.dir_targets_7 = new DenseVector<unsigned long>(dir_targets_7_size);
                fringe.dir_index_7 = new DenseVector<unsigned long>(dir_index_7_size);
                fringe.dir_targets_8 = new DenseVector<unsigned long>(dir_targets_8_size);
                fringe.dir_index_8 = new DenseVector<unsigned long>(dir_index_8_size);
                fringe.external_dir_index_1 = new DenseVector<unsigned long>(external_dir_index_1_size);
                fringe.external_dir_index_2 = new DenseVector<unsigned long>(external_dir_index_2_size);
                fringe.external_dir_index_3 = new DenseVector<unsigned long>(external_dir_index_3_size);
                fringe.external_dir_index_4 = new DenseVector<unsigned long>(external_dir_index_4_size);
                fringe.external_dir_index_5 = new DenseVector<unsigned long>(external_dir_index_5_size);
                fringe.external_dir_index_6 = new DenseVector<unsigned long>(external_dir_index_6_size);
                fringe.external_dir_index_7 = new DenseVector<unsigned long>(external_dir_index_7_size);
                fringe.external_dir_index_8 = new DenseVector<unsigned long>(external_dir_index_8_size);
                fringe.external_h_targets = new DenseVector<unsigned long>(external_h_targets_size);
                fringe.external_dir_targets_1 = new DenseVector<unsigned long>(external_dir_targets_1_size);
                fringe.external_dir_targets_2 = new DenseVector<unsigned long>(external_dir_targets_2_size);
                fringe.external_dir_targets_3 = new DenseVector<unsigned long>(external_dir_targets_3_size);
                fringe.external_dir_targets_4 = new DenseVector<unsigned long>(external_dir_targets_4_size);
                fringe.external_dir_targets_5 = new DenseVector<unsigned long>(external_dir_targets_5_size);
                fringe.external_dir_targets_6 = new DenseVector<unsigned long>(external_dir_targets_6_size);
                fringe.external_dir_targets_7 = new DenseVector<unsigned long>(external_dir_targets_7_size);
                fringe.external_dir_targets_8 = new DenseVector<unsigned long>(external_dir_targets_8_size);

                mpi::mpi_recv(fringe.h_targets->elements(), fringe.h_targets->size(), 0, 0);
                mpi::mpi_recv(fringe.h_index->elements(), fringe.h_index->size(), 0, 0);
                mpi::mpi_recv(fringe.external_h_index->elements(), fringe.external_h_index->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_1->elements(), fringe.dir_targets_1->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_1->elements(), fringe.dir_index_1->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_2->elements(), fringe.dir_targets_2->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_2->elements(), fringe.dir_index_2->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_3->elements(), fringe.dir_targets_3->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_3->elements(), fringe.dir_index_3->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_4->elements(), fringe.dir_targets_4->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_4->elements(), fringe.dir_index_4->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_5->elements(), fringe.dir_targets_5->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_5->elements(), fringe.dir_index_5->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_6->elements(), fringe.dir_targets_6->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_6->elements(), fringe.dir_index_6->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_7->elements(), fringe.dir_targets_7->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_7->elements(), fringe.dir_index_7->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_targets_8->elements(), fringe.dir_targets_8->size(), 0, 0);
                mpi::mpi_recv(fringe.dir_index_8->elements(), fringe.dir_index_8->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_1->elements(), fringe.external_dir_index_1->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_2->elements(), fringe.external_dir_index_2->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_3->elements(), fringe.external_dir_index_3->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_4->elements(), fringe.external_dir_index_4->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_5->elements(), fringe.external_dir_index_5->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_6->elements(), fringe.external_dir_index_6->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_7->elements(), fringe.external_dir_index_7->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_index_8->elements(), fringe.external_dir_index_8->size(), 0, 0);
                mpi::mpi_recv(fringe.external_h_targets->elements(), fringe.external_h_targets->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_1->elements(), fringe.external_dir_targets_1->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_2->elements(), fringe.external_dir_targets_2->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_3->elements(), fringe.external_dir_targets_3->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_4->elements(), fringe.external_dir_targets_4->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_5->elements(), fringe.external_dir_targets_5->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_6->elements(), fringe.external_dir_targets_6->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_7->elements(), fringe.external_dir_targets_7->size(), 0, 0);
                mpi::mpi_recv(fringe.external_dir_targets_8->elements(), fringe.external_dir_targets_8->size(), 0, 0);
            }

            void _send_full_sync(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                mpi::mpi_send(data.h->elements(), data.h->size(), target, _myid);

                mpi::mpi_send(data.f_temp_1->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_2->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_3->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_4->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_5->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_6->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_7->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_8->elements(), data.h->size(), target, _myid);
            }

            void _recv_full_sync(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                mpi::mpi_recv(data.h->elements(), data.h->size(), target, target);

                mpi::mpi_recv(data.f_temp_1->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_2->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_3->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_4->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_5->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_6->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_7->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_8->elements(), data.h->size(), target, target);
            }

            void _send_master_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long offset(info.offset);
                unsigned long f1_offset((*fringe.dir_index_1)[0]);
                unsigned long f1_size((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset);
                if (f1_size > 0) mpi::mpi_send(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, _myid);
                unsigned long f2_offset((*fringe.dir_index_2)[0]);
                unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0)mpi::mpi_send(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, _myid);
                unsigned long f3_offset((*fringe.dir_index_3)[0]);
                unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_send(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, _myid);
                unsigned long f4_offset((*fringe.dir_index_4)[0]);
                unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_send(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, _myid);
                unsigned long f5_offset((*fringe.dir_index_5)[0]);
                unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_send(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, _myid);
                unsigned long f6_offset((*fringe.dir_index_6)[0]);
                unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_send(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, _myid);
                unsigned long f7_offset((*fringe.dir_index_7)[0]);
                unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_send(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, _myid);
                unsigned long f8_offset((*fringe.dir_index_8)[0]);
                unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_send(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, _myid);

                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_send(data.h->elements() + h_offset - offset, h_size, target, _myid);
                }

            }

            void _recv_slave_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long offset(info.offset);
                unsigned long f1_offset((*fringe.dir_index_1)[0]);
                unsigned long f1_size((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset);
                if (f1_size > 0) mpi::mpi_recv(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, target);
                unsigned long f2_offset((*fringe.dir_index_2)[0]);
                unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0) mpi::mpi_recv(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, target);
                unsigned long f3_offset((*fringe.dir_index_3)[0]);
                unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_recv(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, target);
                unsigned long f4_offset((*fringe.dir_index_4)[0]);
                unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_recv(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, target);
                unsigned long f5_offset((*fringe.dir_index_5)[0]);
                unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_recv(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, target);
                unsigned long f6_offset((*fringe.dir_index_6)[0]);
                unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_recv(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, target);
                unsigned long f7_offset((*fringe.dir_index_7)[0]);
                unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_recv(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, target);
                unsigned long f8_offset((*fringe.dir_index_8)[0]);
                unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_recv(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, target);

                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_recv(data.h->elements() + h_offset - offset, h_size, target, target);
                }
            }

            void _send_slave_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long offset(info.offset);
                unsigned long f1_offset((*fringe.external_dir_index_1)[0]);
                unsigned long f1_size((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset);
                if (f1_size > 0) mpi::mpi_send(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, _myid);
                unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0)mpi::mpi_send(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, _myid);
                unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_send(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, _myid);
                unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_send(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, _myid);
                unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_send(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, _myid);
                unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_send(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, _myid);
                unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_send(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, _myid);
                unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_send(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, _myid);

                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_send(data.h->elements() + h_offset - offset, h_size, target, _myid);
                }
            }

            void _recv_master_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                unsigned long offset(info.offset);
                unsigned long f1_offset((*fringe.external_dir_index_1)[0]);
                unsigned long f1_size((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset);
                if (f1_size > 0) mpi::mpi_recv(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, target);
                unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0) mpi::mpi_recv(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, target);
                unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_recv(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, target);
                unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_recv(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, target);
                unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_recv(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, target);
                unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_recv(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, target);
                unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_recv(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, target);
                unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_recv(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, target);

                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_recv(data.h->elements() + h_offset - offset, h_size, target, target);
                }
            }

            void _circle_sync(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                std::vector<MPI::Request> requests;

                unsigned long offset(info.offset);

                {
                    unsigned long source_1((*fringe.external_dir_targets_1)[0] + 1);
                    unsigned long f1_offset((*fringe.external_dir_index_1)[0]);
                    unsigned long f1_size((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset);
                    if (f1_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_1->elements() + f1_offset - offset, f1_size, mpi::MPIType<DataType_>::value(), source_1, source_1));

                    unsigned long source_2((*fringe.external_dir_targets_2)[0] + 1);
                    unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                    unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                    if (f2_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_2->elements() + f2_offset - offset, f2_size, mpi::MPIType<DataType_>::value(), source_2, source_2));

                    unsigned long source_3((*fringe.external_dir_targets_3)[0] + 1);
                    unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                    unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                    if (f3_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_3->elements() + f3_offset - offset, f3_size, mpi::MPIType<DataType_>::value(), source_3, source_3));

                    unsigned long source_4((*fringe.external_dir_targets_4)[0] + 1);
                    unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                    unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                    if (f4_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_4->elements() + f4_offset - offset, f4_size, mpi::MPIType<DataType_>::value(), source_4, source_4));

                    unsigned long source_5((*fringe.external_dir_targets_5)[0] + 1);
                    unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                    unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                    if (f5_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_5->elements() + f5_offset - offset, f5_size, mpi::MPIType<DataType_>::value(), source_5, source_5));

                    unsigned long source_6((*fringe.external_dir_targets_6)[0] + 1);
                    unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                    unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                    if (f6_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_6->elements() + f6_offset - offset, f6_size, mpi::MPIType<DataType_>::value(), source_6, source_6));

                    unsigned long source_7((*fringe.external_dir_targets_7)[0] + 1);
                    unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                    unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                    if (f7_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_7->elements() + f7_offset - offset, f7_size, mpi::MPIType<DataType_>::value(), source_7, source_7));

                    unsigned long source_8((*fringe.external_dir_targets_8)[0] + 1);
                    unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                    unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                    if (f8_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.f_temp_8->elements() + f8_offset - offset, f8_size, mpi::MPIType<DataType_>::value(), source_8, source_8));

                    for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                    {
                        unsigned long h_source((*fringe.h_targets)[i] + 1);
                        unsigned long h_offset((*fringe.h_index)[i * 2]);
                        unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                        if (h_size > 0) requests.push_back(MPI::COMM_WORLD.Irecv(data.h->elements() + h_offset - offset, h_size, mpi::MPIType<DataType_>::value(), h_source, h_source));
                    }
                }

                {
                    unsigned long target_1((*fringe.dir_targets_1)[0] + 1);
                    unsigned long f1_offset((*fringe.dir_index_1)[0]);
                    unsigned long f1_size((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset);
                    if (f1_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_1->elements() + f1_offset - offset, f1_size, mpi::MPIType<DataType_>::value(), target_1, _myid));

                    unsigned long target_2((*fringe.dir_targets_2)[0] + 1);
                    unsigned long f2_offset((*fringe.dir_index_2)[0]);
                    unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                    if (f2_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_2->elements() + f2_offset - offset, f2_size, mpi::MPIType<DataType_>::value(), target_2, _myid));

                    unsigned long target_3((*fringe.dir_targets_3)[0] + 1);
                    unsigned long f3_offset((*fringe.dir_index_3)[0]);
                    unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                    if (f3_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_3->elements() + f3_offset - offset, f3_size, mpi::MPIType<DataType_>::value(), target_3, _myid));

                    unsigned long target_4((*fringe.dir_targets_4)[0] + 1);
                    unsigned long f4_offset((*fringe.dir_index_4)[0]);
                    unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                    if (f4_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_4->elements() + f4_offset - offset, f4_size, mpi::MPIType<DataType_>::value(), target_4, _myid));

                    unsigned long target_5((*fringe.dir_targets_5)[0] + 1);
                    unsigned long f5_offset((*fringe.dir_index_5)[0]);
                    unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                    if (f5_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_5->elements() + f5_offset - offset, f5_size, mpi::MPIType<DataType_>::value(), target_5, _myid));

                    unsigned long target_6((*fringe.dir_targets_6)[0] + 1);
                    unsigned long f6_offset((*fringe.dir_index_6)[0]);
                    unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                    if (f6_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_6->elements() + f6_offset - offset, f6_size, mpi::MPIType<DataType_>::value(), target_6, _myid));

                    unsigned long target_7((*fringe.dir_targets_7)[0] + 1);
                    unsigned long f7_offset((*fringe.dir_index_7)[0]);
                    unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                    if (f7_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_7->elements() + f7_offset - offset, f7_size, mpi::MPIType<DataType_>::value(), target_7, _myid));

                    unsigned long target_8((*fringe.dir_targets_8)[0] + 1);
                    unsigned long f8_offset((*fringe.dir_index_8)[0]);
                    unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                    if (f8_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.f_temp_8->elements() + f8_offset - offset, f8_size, mpi::MPIType<DataType_>::value(), target_8, _myid));

                    for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                    {
                        unsigned long h_target((*fringe.external_h_targets)[i] + 1);
                        unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                        unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                        if (h_size > 0) requests.push_back(MPI::COMM_WORLD.Isend(data.h->elements() + h_offset - offset, h_size, mpi::MPIType<DataType_>::value(), h_target, _myid));
                    }
                }

                MPI::Request::Waitall(requests.size(), &requests[0]);
                requests.clear();
            }
    };
}
#endif
