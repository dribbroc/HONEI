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
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/grid_partitioner.hh>
#include <honei/util/time_stamp.hh>
#include <honei/lbm/scenario_collection.hh>

#include <iostream>
#include <vector>
#include <list>

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
                std::cout<<"Ring LBM Solver with " << _numprocs << " nodes:" << std::endl << std::endl;
                unsigned long timesteps(100);
                Grid<D2Q9, DataType_> grid_global;
                ScenarioCollection::get_scenario(0, gridsize, gridsize, grid_global);
                std::cout << "Solving: " << grid_global.long_description << std::endl;
                std::cout<<"Gridsize: "<<grid_global.h->rows()<<" x "<<grid_global.h->columns()<<std::endl;

                PackedGridData<D2Q9, DataType_>  data_global;
                PackedGridInfo<D2Q9> info_global;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_global, info_global, data_global);
                std::vector<PackedGridInfo<D2Q9> > info_list;
                std::vector<PackedGridData<D2Q9, DataType_> > data_list;
                std::vector<PackedGridFringe<D2Q9> > fringe_list;
                GridPartitioner<D2Q9, DataType_>::decompose(_numprocs, info_global, data_global, info_list, data_list, fringe_list);

                PackedGridData<D2Q9, DataType_> data_lokal(data_list.at(0));
                PackedGridInfo<D2Q9> info_lokal(info_list.at(0));

                mpi::mpi_bcast(&timesteps, 1, 0);
                mpi::mpi_bcast(&grid_global.d_x, 1, 0);
                mpi::mpi_bcast(&grid_global.d_y, 1, 0);
                mpi::mpi_bcast(&grid_global.d_t, 1, 0);
                mpi::mpi_bcast(&grid_global.tau, 1, 0);

                std::list<unsigned long> ul_buffer;
                std::vector<MPI_Request> requests;
                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _send_info(target, info_list.at(target), ul_buffer, requests);
                    _send_data(target, data_list.at(target), ul_buffer, requests);
                    _send_fringe(target, fringe_list.at(target), ul_buffer, requests);
                }
                MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
                requests.clear();
                ul_buffer.clear();


                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info_lokal, &data_lokal, grid_global.d_x, grid_global.d_y, grid_global.d_t, grid_global.tau);

                solver.do_preprocessing();

                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _recv_slave_sync(target, info_list.at(target), data_list.at(target), fringe_list.at(target));
                }

                GridPartitioner<D2Q9, DataType_>::synch(info_global, data_global, info_list, data_list, fringe_list);

                for (signed long target(1) ; target < _numprocs ; ++target)
                {
                    _send_slave_sync(target, info_list.at(target), data_list.at(target), fringe_list.at(target));
                }

                //GridPartitioner<D2Q9, DataType_>::compose(info_global, data_global, info_list, data_list);
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_global, info_global, data_global);
                //PostProcessing<output_types::GNUPLOT>::value(*grid_global.h, 1, grid_global.h->columns(), grid_global.h->rows(), 101);
                // Preproc finished
                // start timesteps
                TimeStamp at, bt;
                at.take();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver.solve();

                    MPI_File fh;
                    std::string filename("h_"+stringify(i)+"_out.dat");
                    MPI_File_open(MPI_COMM_WORLD, (char*)filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
                    MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                    //MPI_File_set_size(fh, 50*50*sizeof(float));
                    MPI_Request request;
                    MPI_File_iwrite_at(fh, info_lokal.offset + (*info_lokal.limits)[0], data_lokal.h->elements() + (*info_lokal.limits)[0], (*info_lokal.limits)[info_lokal.limits->size() - 1] - (*info_lokal.limits)[0], mpi::MPIType<DataType_>::value(), &request);

                    _circle_sync(info_lokal, data_lokal, fringe_list.at(0));

                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                    MPI_File_close(&fh);
                }
                //MPI_Barrier(MPI_COMM_WORLD);
                bt.take();
                std::cout<<"Timesteps: " << timesteps << " TOE: "<<bt.total() - at.total()<<std::endl;
                std::cout<<"MLUPS: "<< (double(grid_global.h->rows()) * double(grid_global.h->columns()) * double(timesteps)) / (1e6 * (bt.total() - at.total())) <<std::endl;

                /*std::cout<<*data_global.h;
                MPI_File fh;
                std::string filename("h_"+stringify(99)+"_out.dat");
                MPI_File_open(MPI_COMM_SELF, (char*)filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                MPI_File_read_at(fh, 0, data_global.h->elements(), (*info_global.limits)[info_global.limits->size() - 1] - (*info_global.limits)[0], mpi::MPIType<DataType_>::value(), MPI_STATUSES_IGNORE);
                MPI_File_close(&fh);
                std::cout<<*data_global.h;*/

                // generate output
                /*for (signed long target(1) ; target < _numprocs ; ++target)
                  {
                  _recv_full_sync(target, data_list.at(target));
                }
                GridPartitioner<D2Q9, DataType_>::synch(info_global, data_global, info_list, data_list, fringe_list);
                GridPartitioner<D2Q9, DataType_>::compose(info_global, data_global, info_list, data_list);*/
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_global, info_global, data_global);
                //std::cout<<*grid_global.h;
                //PostProcessing<output_types::GNUPLOT>::value(*grid_global.h, 1, grid_global.h->columns(), grid_global.h->rows(),

                /*Grid<D2Q9, DataType_> grid_ref;
                ScenarioCollection::get_scenario(0, gridsize, gridsize, grid_ref);

                PackedGridData<D2Q9, DataType_>  data_ref;
                PackedGridInfo<D2Q9> info_ref;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_ref, info_ref, data_ref);
                //SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info_ref, &data_ref, grid_ref.d_x, grid_ref.d_y, grid_ref.d_t, grid_ref.tau);
                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver_ref(&info_ref, &data_ref, grid_ref.d_x, grid_ref.d_y, grid_ref.d_t, grid_ref.tau);
                solver_ref.do_preprocessing();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver_ref.solve();
                }
                data_ref.h->lock(lm_read_only);
                for (unsigned long i(0) ; i < data_global.h->size() ; ++i)
                {
                    if (fabs((*data_global.h)[i] - (*data_ref.h)[i]) > 0.0001)
                        std::cout<<(*data_global.h)[i]<<" "<<(*data_ref.h)[i]<<std::endl;
                }
                data_ref.h->unlock(lm_read_only);*/
            }

            void _slave()
            {
                PackedGridData<D2Q9, DataType_> data;
                PackedGridInfo<D2Q9> info;
                PackedGridFringe<D2Q9> fringe;
                unsigned long timesteps;
                DataType_ d_x, d_y, d_t, tau;

                mpi::mpi_bcast(&timesteps, 1, 0);
                mpi::mpi_bcast(&d_x, 1, 0);
                mpi::mpi_bcast(&d_y, 1, 0);
                mpi::mpi_bcast(&d_t, 1, 0);
                mpi::mpi_bcast(&tau, 1, 0);

                _recv_info(info);
                _recv_data(data);
                _recv_fringe(fringe);

                //SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, d_x, d_y, d_t, tau);
                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, d_x, d_y, d_t, tau);

                solver.do_preprocessing();

                _send_master_sync(0, info, data, fringe);
                _recv_master_sync(0, info, data, fringe);

                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver.solve();

                    MPI_File fh;
                    std::string filename("h_"+stringify(i)+"_out.dat");
                    MPI_File_open(MPI_COMM_WORLD, (char*)filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
                    MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                    //MPI_File_set_size(fh, 50*50*sizeof(float));
                    MPI_Request request;
                    MPI_File_iwrite_at(fh, info.offset + (*info.limits)[0], data.h->elements() + (*info.limits)[0], (*info.limits)[info.limits->size() - 1] - (*info.limits)[0], mpi::MPIType<DataType_>::value(), &request);

                    _circle_sync(info, data, fringe);

                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                    MPI_File_close(&fh);
                }
                //MPI_Barrier(MPI_COMM_WORLD);
                //_send_full_sync(0, data);
            }

            void _send_info(unsigned long target, PackedGridInfo<D2Q9> & info, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(info.offset);
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.limits->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(info.dir_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));

                requests.push_back(mpi::mpi_isend(info.limits->elements(), info.limits->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.types->elements(), info.types->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_1->elements(), info.dir_index_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_1->elements(), info.dir_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_2->elements(), info.dir_index_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_2->elements(), info.dir_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_3->elements(), info.dir_index_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_3->elements(), info.dir_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_4->elements(), info.dir_index_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_4->elements(), info.dir_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_5->elements(), info.dir_index_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_5->elements(), info.dir_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_6->elements(), info.dir_index_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_6->elements(), info.dir_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_7->elements(), info.dir_index_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_7->elements(), info.dir_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_index_8->elements(), info.dir_index_8->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(info.dir_8->elements(), info.dir_8->size(), target, _myid));
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

            void _send_data(unsigned long target, PackedGridData<D2Q9, DataType_> & data, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(data.h->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(data.distribution_x->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));

                requests.push_back(mpi::mpi_isend(data.h->elements(), data.h->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(data.u->elements(), data.u->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(data.v->elements(), data.v->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(data.b->elements(), data.b->size(), target, _myid));
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

            void _send_fringe(unsigned long target, PackedGridFringe<D2Q9> & fringe, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(fringe.h_targets->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.h_index->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_h_index->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_targets_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_h_targets->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));
                ul_buffer.push_back(fringe.external_dir_targets_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _myid));

                requests.push_back(mpi::mpi_isend(fringe.h_targets->elements(), fringe.h_targets->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.h_index->elements(), fringe.h_index->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_h_index->elements(), fringe.external_h_index->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_1->elements(), fringe.dir_targets_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_1->elements(), fringe.dir_index_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_2->elements(), fringe.dir_targets_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_2->elements(), fringe.dir_index_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_3->elements(), fringe.dir_targets_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_3->elements(), fringe.dir_index_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_4->elements(), fringe.dir_targets_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_4->elements(), fringe.dir_index_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_5->elements(), fringe.dir_targets_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_5->elements(), fringe.dir_index_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_6->elements(), fringe.dir_targets_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_6->elements(), fringe.dir_index_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_7->elements(), fringe.dir_targets_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_7->elements(), fringe.dir_index_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_8->elements(), fringe.dir_targets_8->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_8->elements(), fringe.dir_index_8->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_1->elements(), fringe.external_dir_index_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_2->elements(), fringe.external_dir_index_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_3->elements(), fringe.external_dir_index_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_4->elements(), fringe.external_dir_index_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_5->elements(), fringe.external_dir_index_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_6->elements(), fringe.external_dir_index_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_7->elements(), fringe.external_dir_index_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_8->elements(), fringe.external_dir_index_8->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_h_targets->elements(), fringe.external_h_targets->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_1->elements(), fringe.external_dir_targets_1->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_2->elements(), fringe.external_dir_targets_2->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_3->elements(), fringe.external_dir_targets_3->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_4->elements(), fringe.external_dir_targets_4->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_5->elements(), fringe.external_dir_targets_5->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_6->elements(), fringe.external_dir_targets_6->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_7->elements(), fringe.external_dir_targets_7->size(), target, _myid));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_8->elements(), fringe.external_dir_targets_8->size(), target, _myid));
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
                data.h->lock(lm_read_only);
                data.f_temp_1->lock(lm_read_only);
                data.f_temp_2->lock(lm_read_only);
                data.f_temp_3->lock(lm_read_only);
                data.f_temp_4->lock(lm_read_only);
                data.f_temp_5->lock(lm_read_only);
                data.f_temp_6->lock(lm_read_only);
                data.f_temp_7->lock(lm_read_only);
                data.f_temp_8->lock(lm_read_only);

                mpi::mpi_send(data.h->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_1->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_2->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_3->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_4->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_5->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_6->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_7->elements(), data.h->size(), target, _myid);
                mpi::mpi_send(data.f_temp_8->elements(), data.h->size(), target, _myid);

                data.h->unlock(lm_read_only);
                data.f_temp_1->unlock(lm_read_only);
                data.f_temp_2->unlock(lm_read_only);
                data.f_temp_3->unlock(lm_read_only);
                data.f_temp_4->unlock(lm_read_only);
                data.f_temp_5->unlock(lm_read_only);
                data.f_temp_6->unlock(lm_read_only);
                data.f_temp_7->unlock(lm_read_only);
                data.f_temp_8->unlock(lm_read_only);
            }

            void _recv_full_sync(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                data.h->lock(lm_read_and_write);
                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

                mpi::mpi_recv(data.h->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_1->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_2->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_3->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_4->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_5->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_6->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_7->elements(), data.h->size(), target, target);
                mpi::mpi_recv(data.f_temp_8->elements(), data.h->size(), target, target);

                data.h->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }

            void _send_master_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                data.h->lock(lm_read_only);
                data.f_temp_1->lock(lm_read_only);
                data.f_temp_2->lock(lm_read_only);
                data.f_temp_3->lock(lm_read_only);
                data.f_temp_4->lock(lm_read_only);
                data.f_temp_5->lock(lm_read_only);
                data.f_temp_6->lock(lm_read_only);
                data.f_temp_7->lock(lm_read_only);
                data.f_temp_8->lock(lm_read_only);

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

                data.h->unlock(lm_read_only);
                data.f_temp_1->unlock(lm_read_only);
                data.f_temp_2->unlock(lm_read_only);
                data.f_temp_3->unlock(lm_read_only);
                data.f_temp_4->unlock(lm_read_only);
                data.f_temp_5->unlock(lm_read_only);
                data.f_temp_6->unlock(lm_read_only);
                data.f_temp_7->unlock(lm_read_only);
                data.f_temp_8->unlock(lm_read_only);
            }

            void _recv_slave_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                data.h->lock(lm_read_and_write);
                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

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

                data.h->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }

            void _send_slave_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                data.h->lock(lm_read_only);
                data.f_temp_1->lock(lm_read_only);
                data.f_temp_2->lock(lm_read_only);
                data.f_temp_3->lock(lm_read_only);
                data.f_temp_4->lock(lm_read_only);
                data.f_temp_5->lock(lm_read_only);
                data.f_temp_6->lock(lm_read_only);
                data.f_temp_7->lock(lm_read_only);
                data.f_temp_8->lock(lm_read_only);

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

                data.h->unlock(lm_read_only);
                data.f_temp_1->unlock(lm_read_only);
                data.f_temp_2->unlock(lm_read_only);
                data.f_temp_3->unlock(lm_read_only);
                data.f_temp_4->unlock(lm_read_only);
                data.f_temp_5->unlock(lm_read_only);
                data.f_temp_6->unlock(lm_read_only);
                data.f_temp_7->unlock(lm_read_only);
                data.f_temp_8->unlock(lm_read_only);
            }

            void _recv_master_sync(unsigned long target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                data.h->lock(lm_read_and_write);
                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

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

                data.h->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }

            void _circle_sync(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
            {
                /// \todo global buffer for requests and in/out data

                data.h->lock(lm_read_and_write);
                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

                std::vector<MPI_Request> requests;
                unsigned long offset(info.offset);

                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    unsigned long h_source((*fringe.h_targets)[i]);
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) requests.push_back(mpi::mpi_irecv(data.h->elements() + h_offset - offset, h_size, h_source, h_source));
                }
                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    unsigned long h_target((*fringe.external_h_targets)[i]);
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) requests.push_back(mpi::mpi_isend(data.h->elements() + h_offset - offset, h_size, h_target, _myid));
                }

                unsigned long source_1_recv((*fringe.external_dir_targets_1)[0]);
                unsigned long f1_offset_recv((*fringe.external_dir_index_1)[0]);
                unsigned long f1_size_recv((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset_recv);

                unsigned long source_2_recv((*fringe.external_dir_targets_2)[0]);
                unsigned long f2_offset_recv((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size_recv((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset_recv);

                unsigned long source_3_recv((*fringe.external_dir_targets_3)[0]);
                unsigned long f3_offset_recv((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size_recv((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset_recv);

                unsigned long source_4_recv((*fringe.external_dir_targets_4)[0]);
                unsigned long f4_offset_recv((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size_recv((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset_recv);

                unsigned long source_5_recv((*fringe.external_dir_targets_5)[0]);
                unsigned long f5_offset_recv((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size_recv((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset_recv);

                unsigned long source_6_recv((*fringe.external_dir_targets_6)[0]);
                unsigned long f6_offset_recv((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size_recv((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset_recv);

                unsigned long source_7_recv((*fringe.external_dir_targets_7)[0]);
                unsigned long f7_offset_recv((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size_recv((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset_recv);

                unsigned long source_8_recv((*fringe.external_dir_targets_8)[0]);
                unsigned long f8_offset_recv((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size_recv((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset_recv);

                unsigned long up_size_recv(f2_size_recv + f3_size_recv + f4_size_recv + f5_size_recv);
                DataType_ up_buffer_recv[up_size_recv];
                unsigned long source_up_recv(source_2_recv);
                source_up_recv = std::max(source_up_recv, source_3_recv);
                source_up_recv = std::max(source_up_recv, source_4_recv);
                source_up_recv = std::max(source_up_recv, source_5_recv);
                unsigned long down_size_recv(f1_size_recv + f6_size_recv + f7_size_recv + f8_size_recv);
                DataType_ down_buffer_recv[down_size_recv];
                unsigned long source_down_recv(source_1_recv);
                source_down_recv = std::max(source_down_recv, source_6_recv);
                source_down_recv = std::max(source_down_recv, source_7_recv);
                source_down_recv = std::max(source_down_recv, source_8_recv);

                if (up_size_recv > 0) requests.push_back(mpi::mpi_irecv(up_buffer_recv, up_size_recv, source_up_recv, source_up_recv));
                if (down_size_recv > 0) requests.push_back(mpi::mpi_irecv(down_buffer_recv, down_size_recv, source_down_recv, source_down_recv));


                unsigned long target_1_send((*fringe.dir_targets_1)[0]);
                unsigned long f1_offset_send((*fringe.dir_index_1)[0]);
                unsigned long f1_size_send((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset_send);

                unsigned long target_2_send((*fringe.dir_targets_2)[0]);
                unsigned long f2_offset_send((*fringe.dir_index_2)[0]);
                unsigned long f2_size_send((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset_send);

                unsigned long target_3_send((*fringe.dir_targets_3)[0]);
                unsigned long f3_offset_send((*fringe.dir_index_3)[0]);
                unsigned long f3_size_send((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset_send);

                unsigned long target_4_send((*fringe.dir_targets_4)[0]);
                unsigned long f4_offset_send((*fringe.dir_index_4)[0]);
                unsigned long f4_size_send((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset_send);

                unsigned long target_5_send((*fringe.dir_targets_5)[0]);
                unsigned long f5_offset_send((*fringe.dir_index_5)[0]);
                unsigned long f5_size_send((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset_send);

                unsigned long target_6_send((*fringe.dir_targets_6)[0]);
                unsigned long f6_offset_send((*fringe.dir_index_6)[0]);
                unsigned long f6_size_send((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset_send);

                unsigned long target_7_send((*fringe.dir_targets_7)[0]);
                unsigned long f7_offset_send((*fringe.dir_index_7)[0]);
                unsigned long f7_size_send((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset_send);

                unsigned long target_8_send((*fringe.dir_targets_8)[0]);
                unsigned long f8_offset_send((*fringe.dir_index_8)[0]);
                unsigned long f8_size_send((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset_send);


                unsigned long up_size_send(f2_size_send + f3_size_send + f4_size_send + f5_size_send);
                DataType_ up_buffer_send[up_size_send];
                unsigned long target_up_send(target_2_send);
                target_up_send = std::max(target_up_send, target_3_send);
                target_up_send = std::max(target_up_send, target_4_send);
                target_up_send = std::max(target_up_send, target_5_send);
                unsigned long down_size_send(f1_size_send + f6_size_send + f7_size_send + f8_size_send);
                DataType_ down_buffer_send[down_size_send];
                unsigned long target_down_send(target_1_send);
                target_down_send = std::max(target_down_send, target_6_send);
                target_down_send = std::max(target_down_send, target_7_send);
                target_down_send = std::max(target_down_send, target_8_send);
                unsigned long temp_size(0);

                TypeTraits<DataType_>::copy(data.f_temp_2->elements() + f2_offset_send - offset, up_buffer_send + temp_size, f2_size_send);
                temp_size += f2_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_3->elements() + f3_offset_send - offset, up_buffer_send + temp_size, f3_size_send);
                temp_size += f3_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_4->elements() + f4_offset_send - offset, up_buffer_send + temp_size, f4_size_send);
                temp_size += f4_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_5->elements() + f5_offset_send - offset, up_buffer_send + temp_size, f5_size_send);

                temp_size = 0;
                TypeTraits<DataType_>::copy(data.f_temp_1->elements() + f1_offset_send - offset, down_buffer_send + temp_size, f1_size_send);
                temp_size += f1_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_6->elements() + f6_offset_send - offset, down_buffer_send + temp_size, f6_size_send);
                temp_size += f6_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_7->elements() + f7_offset_send - offset, down_buffer_send + temp_size, f7_size_send);
                temp_size += f7_size_send;
                TypeTraits<DataType_>::copy(data.f_temp_8->elements() + f8_offset_send - offset, down_buffer_send + temp_size, f8_size_send);

                if (up_size_send > 0) requests.push_back(mpi::mpi_isend(up_buffer_send, up_size_send, target_up_send, _myid));
                if (down_size_send > 0) requests.push_back(mpi::mpi_isend(down_buffer_send, down_size_send, target_down_send, _myid));


                MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

                temp_size = 0;
                TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_2->elements() + f2_offset_recv - offset, f2_size_recv);
                temp_size += f2_size_recv;
                TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_3->elements() + f3_offset_recv - offset, f3_size_recv);
                temp_size += f3_size_recv;
                TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_4->elements() + f4_offset_recv - offset, f4_size_recv);
                temp_size += f4_size_recv;
                TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_5->elements() + f5_offset_recv - offset, f5_size_recv);

                temp_size = 0;
                TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_1->elements() + f1_offset_recv - offset, f1_size_recv);
                temp_size += f1_size_recv;
                TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_6->elements() + f6_offset_recv - offset, f6_size_recv);
                temp_size += f6_size_recv;
                TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_7->elements() + f7_offset_recv - offset, f7_size_recv);
                temp_size += f7_size_recv;
                TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_8->elements() + f8_offset_recv - offset, f8_size_recv);

                data.h->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };
}
#endif
