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
#include <honei/math/vector_io.hh>
#include <honei/util/string_tokenizer.hh>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
//#include <netdb.h>

namespace honei
{
    template <typename DataType_>
    class MPIRingSolver
    {
        private:
            int _numprocs;
            int _nodes;
            int _myid;
            int _mycartid;
            int _masterid;
            int _mycartpos;
            MPI_Comm _comm_cart;
            std::string _base_file_name;
            bool _file_output;
            int _gpu_device;
            tags::TagValue _solver_tag_value;
            std::vector<std::string> _backends;
            std::vector<double> _fractions;
            double _sync_time_up;
            double _sync_time_down;
            std::string _device_name;
            unsigned long _scenario;
            double _sync_threshold;
            bool _recently_synched;

        public:
            MPIRingSolver(int argc, char **argv)
            {
                mpi::mpi_init(&argc, &argv);
                mpi::mpi_comm_size(&_numprocs);
                mpi::mpi_comm_rank(&_myid);
                _gpu_device = 0;
                _sync_time_up = 0;
                _sync_time_down = 0;
                _recently_synched = false;
                if (argc != 5)
                {
                    if(_myid == 0) std::cout<<"Usage: honei-mpi-ring grid_x grid_y timesteps config_file_name"<<std::endl;
                    mpi::mpi_finalize();
                    exit(1);
                }
                unsigned long gridsize_x(atoi(argv[1]));
                unsigned long gridsize_y(atoi(argv[2]));
                unsigned long timesteps(atoi(argv[3]));
                std::string config_file_name(argv[4]);
                _file_output = false;

                // create new communicator
                int dims[1];
                dims[0] = _numprocs;
                int periods[1];
                periods[0] = false;
                MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, true, &_comm_cart);
                mpi::mpi_comm_rank(&_mycartid, _comm_cart);
                int position(0);
                MPI_Cart_rank(_comm_cart, &position, &_masterid);
                MPI_Cart_coords(_comm_cart, _mycartid, 1, &_mycartpos);
                /*std::cout<<_myid<<" "<<_mycartid<<std::endl;
                int test;
                int input(4);
                MPI_Cart_rank(_comm_cart, &input, &test);
                std::cout<<test<<std::endl;;
                MPI_Cart_shift(_comm_cart, 0, 1, &input, &test);
                std::cout<<input<<" -> "<<test<<std::endl;*/

                //char hostname [256];
                //gethostname(hostname, 255);
                //std::cout<<"this is process " << _mycartid << " on machine "<<hostname<<std::endl;

                //read in configuration file
                _read_config(config_file_name, _scenario, _backends, _fractions, _file_output, _base_file_name);
                // \TODO default werte setzen falls kein config eintrag?

                _nodes = _numprocs / (_fractions.size());


                if (_mycartid == _masterid)
                {
                    _master(gridsize_x, gridsize_y, timesteps);
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
            void _master(unsigned long gridsize_x, unsigned long gridsize_y, unsigned long timesteps)
            {
                std::cout<<"Ring LBM Solver with " << _numprocs << " procs on " << _nodes << " nodes (" << _numprocs / _nodes << " jobs per node)" << std::endl;
                if (_file_output)
                    std::cout<<"with file output activated"<<std::endl<<std::endl;
                else
                    std::cout<<std::endl;
                Grid<D2Q9, DataType_> grid_global;
                ScenarioCollection::get_scenario(_scenario, gridsize_x, gridsize_y, grid_global);
                std::cout << "Solving: " << grid_global.long_description << std::endl;
                std::cout<<"Gridsize: "<<grid_global.h->rows()<<" x "<<grid_global.h->columns()<<std::endl;
                std::cout<<"Timesteps: "<<timesteps<<std::endl;

                PackedGridData<D2Q9, DataType_>  data_global;
                PackedGridInfo<D2Q9> info_global;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_global, info_global, data_global, false);
                std::vector<PackedGridInfo<D2Q9> > info_list;
                std::vector<PackedGridData<D2Q9, DataType_> > data_list;
                std::vector<PackedGridFringe<D2Q9> > fringe_list;

                //decompose patches
                {
                    if (_numprocs % (_fractions.size()) != 0)
                        throw InternalError("numprocs / _fractions missmatch!");

                    std::vector<unsigned long> patch_sizes;
                    unsigned long size_per_node(data_global.u->size() / _nodes);

                    for (long node(0) ; node < _nodes ; ++node)
                    {
                        for (unsigned long i(0) ; i < _fractions.size() ; ++i)
                            patch_sizes.push_back(size_per_node * _fractions.at(i));
                        unsigned long whole_size(0);
                        for (unsigned long i(node * _fractions.size()) ; i < patch_sizes.size() ; ++i)
                            whole_size += patch_sizes.at(i);
                        patch_sizes.at(node * _fractions.size()) += size_per_node - whole_size;
                    }
                    patch_sizes.front() += data_global.u->size() % _nodes;
                    GridPartitioner<D2Q9, DataType_>::decompose_intern(patch_sizes, info_global, data_global, info_list, data_list, fringe_list, false);
                }

                PackedGridData<D2Q9, DataType_> data_lokal(data_list.at(0));
                PackedGridInfo<D2Q9> info_lokal(info_list.at(0));
                _alloc_lokal_data(data_lokal);

                unsigned long mpi_file_size(data_global.h->size() * sizeof(DataType_));
                mpi::mpi_bcast(&timesteps, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&mpi_file_size, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&grid_global.d_x, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&grid_global.d_y, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&grid_global.d_t, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&grid_global.tau, 1, _masterid, _comm_cart);

                std::list<unsigned long> ul_buffer;
                std::vector<MPI_Request> requests;
                for (int target(1) ; target < _numprocs ; ++target)
                {
                    int rank;
                    MPI_Cart_rank(_comm_cart, &target, &rank);
                    _send_info(rank, info_list.at(target), ul_buffer, requests);
                    _send_data(rank, data_list.at(target), ul_buffer, requests);
                    _send_fringe(rank, fringe_list.at(target), ul_buffer, requests);
                }
                MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
                requests.clear();
                ul_buffer.clear();


                SolverLBMGridBase * solver(NULL);
                _init_solver(solver, _mycartid % _backends.size(), info_lokal, data_lokal, grid_global.d_x, grid_global.d_y, grid_global.d_t, grid_global.tau);

                solver->do_preprocessing();

                /*for (int target(1) ; target < _numprocs ; ++target)
                  {
                  int rank;
                  MPI_Cart_rank(_comm_cart, &target, &rank);
                  _recv_slave_sync(rank, info_list.at(target), data_list.at(target), fringe_list.at(target));
                  }

                  GridPartitioner<D2Q9, DataType_>::synch(info_global, data_global, info_list, data_list, fringe_list);

                  for (int target(1) ; target < _numprocs ; ++target)
                  {
                  int rank;
                  MPI_Cart_rank(_comm_cart, &target, &rank);
                  _send_slave_sync(rank, info_list.at(target), data_list.at(target), fringe_list.at(target));
                  }*/

                MPI_Barrier(MPI_COMM_WORLD);
                _circle_sync(info_lokal, data_lokal, fringe_list.at(0));

                //GridPartitioner<D2Q9, DataType_>::compose(info_global, data_global, info_list, data_list);
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_global, info_global, data_global);
                //PostProcessing<output_types::GNUPLOT>::value(*grid_global.h, 1, grid_global.h->columns(), grid_global.h->rows(), 101);
                // Preproc finished
                // start timesteps
                TimeStamp at, bt;
                at.take();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver->solve();

                    /*MPI_File fh;
                      std::string filename("h_"+stringify(i)+"_out.dat");
                      MPI_File_open(MPI_COMM_WORLD, (char*)filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
                      MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                      MPI_File_set_size(fh, mpi_file_size);
                      MPI_Request request;
                      MPI_File_iwrite_at(fh, info_lokal.offset + (*info_lokal.limits)[0], data_lokal.h->elements() + (*info_lokal.limits)[0], (*info_lokal.limits)[info_lokal.limits->size() - 1] - (*info_lokal.limits)[0], mpi::MPIType<DataType_>::value(), &request);*/

                    _circle_sync(info_lokal, data_lokal, fringe_list.at(0));

                    if (_file_output)
                    {
                        std::string fn("h_"+_base_file_name+"_"+stringify(i)+"_"+stringify(_mycartpos)+".dat");
                        DenseVectorRange<DataType_> range(data_lokal.h->range((*info_lokal.limits)[info_lokal.limits->size() - 1] - (*info_lokal.limits)[0], (*info_lokal.limits)[0]));
                        VectorIO<io_formats::DV>::write_vector(fn, range);
                    }

                    //MPI_Wait(&request, MPI_STATUS_IGNORE);
                    //MPI_File_close(&fh);
                }
                std::cout<<_mycartid << " (" << _device_name << "): up "<<_sync_time_up<<" down "<<_sync_time_down<<std::endl;
                MPI_Barrier(MPI_COMM_WORLD);
                bt.take();
                solver->do_postprocessing();
                std::cout<<"Timesteps: " << timesteps << " TOE: "<<bt.total() - at.total()<<std::endl;
                std::cout<<"MLUPS: "<< (double(grid_global.h->rows()) * double(grid_global.h->columns()) * double(timesteps)) / (1e6 * (bt.total() - at.total())) <<std::endl;



                if (_file_output)
                {
                    for (unsigned long i(0) ; i < timesteps ; ++i)
                    {
                        DenseVector<DataType_> global_h(data_global.h->size());
                        unsigned long global_i(0);
                        for (long j(0) ; j < _numprocs ; ++j)
                        {
                            std::string fn("h_"+_base_file_name+"_"+stringify(i)+"_"+stringify(j)+".dat");
                            DenseVector<DataType_> temp(VectorIO<io_formats::DV>::read_vector(fn, DataType_(0)));
                            for (unsigned long x(0) ; x < temp.size() ; ++x)
                            {
                                global_h[global_i] = temp[x];
                                ++global_i;
                            }
                            remove(fn.c_str());
                        }
                        std::string out_filename("h_"+_base_file_name+"_"+stringify(i)+".dv");
                        VectorIO<io_formats::DV>::write_vector(out_filename, global_h);
                    }
                }

                /*std::cout<<*data_global.h;
                  MPI_File fh;
                  std::string filename("h_"+stringify(99)+"_out.dat");
                  MPI_File_open(MPI_COMM_SELF, (char*)filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                  MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                  MPI_File_read_at(fh, 0, data_global.h->elements(), (*info_global.limits)[info_global.limits->size() - 1] - (*info_global.limits)[0], mpi::MPIType<DataType_>::value(), MPI_STATUSES_IGNORE);
                  MPI_File_close(&fh);
                  std::cout<<*data_global.h;*/

                // generate output
                /*for (int target(1) ; target < _numprocs ; ++target)
                {
                    int rank;
                    MPI_Cart_rank(_comm_cart, &target, &rank);
                    _recv_full_sync(target, data_list.at(target));
                }
                GridPartitioner<D2Q9, DataType_>::synch(info_global, data_global, info_list, data_list, fringe_list);
                GridPartitioner<D2Q9, DataType_>::compose(info_global, data_global, info_list, data_list);*/
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_global, info_global, data_global);
                //std::cout<<*grid_global.h;
                //PostProcessing<output_types::GNUPLOT>::value(*grid_global.h, 1, grid_global.h->columns(), grid_global.h->rows(),

                /*Grid<D2Q9, DataType_> grid_ref;
                ScenarioCollection::get_scenario(0, gridsize_x, gridsize_y, grid_ref);

                PackedGridData<D2Q9, DataType_>  data_ref;
                PackedGridInfo<D2Q9> info_ref;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_ref, info_ref, data_ref);
                SolverLBMGrid<tags::CPU::SSE, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver_ref(&info_ref, &data_ref, grid_ref.d_x, grid_ref.d_y, grid_ref.d_t, grid_ref.tau);
                solver_ref.do_preprocessing();
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver_ref.solve();
                }
                solver_ref.do_postprocessing();
                data_ref.h->lock(lm_read_only);
                for (unsigned long i(0) ; i < data_global.h->size() ; ++i)
                {
                    if (fabs((*data_global.h)[i] - (*data_ref.h)[i]) > 0.0001)
                        std::cout<<(*data_global.h)[i]<<" "<<(*data_ref.h)[i]<<std::endl;
                }
                data_ref.h->unlock(lm_read_only);
                grid_ref.destroy();
                info_ref.destroy();
                data_ref.destroy();
                std::cout<<"Comparison finished sucessfully!"<<std::endl;*/

                delete solver;
                data_lokal.destroy();
                data_list.erase(data_list.begin());
                honei::GridPartitioner<D2Q9, DataType_>::destroy(info_list, data_list, fringe_list);
                grid_global.destroy();
                info_global.destroy();
                data_global.destroy();
            }

            void _slave()
            {
                PackedGridData<D2Q9, DataType_> data;
                PackedGridInfo<D2Q9> info;
                PackedGridFringe<D2Q9> fringe;
                unsigned long timesteps, mpi_file_size;
                DataType_ d_x, d_y, d_t, tau;

                mpi::mpi_bcast(&timesteps, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&mpi_file_size, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&d_x, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&d_y, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&d_t, 1, _masterid, _comm_cart);
                mpi::mpi_bcast(&tau, 1, _masterid, _comm_cart);

                _recv_info(info);
                _recv_data(data);
                _recv_fringe(fringe);

                SolverLBMGridBase * solver(NULL);
                _init_solver(solver, _mycartid % _backends.size(), info, data, d_x, d_y, d_t, tau);

                solver->do_preprocessing();

                //_send_master_sync(_masterid, info, data, fringe);
                //_recv_master_sync(_masterid, info, data, fringe);
                MPI_Barrier(MPI_COMM_WORLD);
                _circle_sync(info, data, fringe);

                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver->solve();

                    /*MPI_File fh;
                      std::string filename("h_"+stringify(i)+"_out.dat");
                      MPI_File_open(MPI_COMM_WORLD, (char*)filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
                      MPI_File_set_view(fh, 0, mpi::MPIType<DataType_>::value(), mpi::MPIType<DataType_>::value(), "native", MPI_INFO_NULL);
                      MPI_File_set_size(fh, mpi_file_size);
                      MPI_Request request;
                      MPI_File_iwrite_at(fh, info.offset + (*info.limits)[0], data.h->elements() + (*info.limits)[0], (*info.limits)[info.limits->size() - 1] - (*info.limits)[0], mpi::MPIType<DataType_>::value(), &request);*/

                    _circle_sync(info, data, fringe);

                    if (_file_output)
                    {
                        std::string fn("h_"+_base_file_name+"_"+stringify(i)+"_"+stringify(_mycartpos)+".dat");
                        DenseVectorRange<DataType_> range(data.h->range((*info.limits)[info.limits->size() - 1] - (*info.limits)[0], (*info.limits)[0]));
                        VectorIO<io_formats::DV>::write_vector(fn, range);
                    }

                    //MPI_Wait(&request, MPI_STATUS_IGNORE);
                    //MPI_File_close(&fh);
                }
                std::cout<<_mycartid << " (" << _device_name << "): up "<<_sync_time_up<<" down "<<_sync_time_down<<std::endl;
                MPI_Barrier(MPI_COMM_WORLD);
                solver->do_postprocessing();
                //_send_full_sync(_masterid, data);

                delete solver;
                fringe.destroy();
                info.destroy();
                data.destroy();
            }

            void _alloc_lokal_data(PackedGridData<D2Q9, DataType_> & data_lokal)
            {
                data_lokal.temp = new DenseVector<DataType_>(data_lokal.h->size());

                data_lokal.f_0 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_1 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_2 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_3 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_4 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_5 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_6 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_7 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_8 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));

                data_lokal.f_eq_0 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_1 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_2 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_3 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_4 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_5 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_6 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_7 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_eq_8 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));

                data_lokal.f_temp_0 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_1 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_2 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_3 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_4 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_5 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_6 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_7 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));
                data_lokal.f_temp_8 = new DenseVector<DataType_>(data_lokal.h->size(), DataType_(0));

                data_lokal.distribution_x = new DenseVector<DataType_>(9ul, DataType_(0));
                data_lokal.distribution_y = new DenseVector<DataType_>(9ul, DataType_(0));
            }

            void _send_info(int target, PackedGridInfo<D2Q9> & info, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(info.offset);
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.limits->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(info.dir_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));

                requests.push_back(mpi::mpi_isend(info.limits->elements(), info.limits->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.types->elements(), info.types->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_1->elements(), info.dir_index_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_1->elements(), info.dir_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_2->elements(), info.dir_index_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_2->elements(), info.dir_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_3->elements(), info.dir_index_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_3->elements(), info.dir_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_4->elements(), info.dir_index_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_4->elements(), info.dir_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_5->elements(), info.dir_index_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_5->elements(), info.dir_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_6->elements(), info.dir_index_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_6->elements(), info.dir_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_7->elements(), info.dir_index_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_7->elements(), info.dir_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_index_8->elements(), info.dir_index_8->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(info.dir_8->elements(), info.dir_8->size(), target, _mycartid, _comm_cart));
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

                mpi::mpi_recv(&info.offset, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&limits_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_8_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_8_size, 1, _masterid, _masterid, _comm_cart);

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

                mpi::mpi_recv(info.limits->elements(), info.limits->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.types->elements(), info.types->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_1->elements(), info.dir_index_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_1->elements(), info.dir_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_2->elements(), info.dir_index_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_2->elements(), info.dir_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_3->elements(), info.dir_index_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_3->elements(), info.dir_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_4->elements(), info.dir_index_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_4->elements(), info.dir_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_5->elements(), info.dir_index_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_5->elements(), info.dir_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_6->elements(), info.dir_index_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_6->elements(), info.dir_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_7->elements(), info.dir_index_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_7->elements(), info.dir_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_index_8->elements(), info.dir_index_8->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(info.dir_8->elements(), info.dir_8->size(), _masterid, _masterid, _comm_cart);
            }

            void _send_data(int target, PackedGridData<D2Q9, DataType_> & data, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(data.h->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                //ul_buffer.push_back(data.distribution_x->size());
                ul_buffer.push_back(9);
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));

                requests.push_back(mpi::mpi_isend(data.h->elements(), data.h->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(data.u->elements(), data.u->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(data.v->elements(), data.v->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(data.b->elements(), data.b->size(), target, _mycartid, _comm_cart));
            }

            void _recv_data(PackedGridData<D2Q9, DataType_> & data)
            {
                unsigned long h_size;
                unsigned long dist_size;

                mpi::mpi_recv(&h_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dist_size, 1, _masterid, _masterid, _comm_cart);

                data.h = new DenseVector<DataType_>(h_size);
                data.u = new DenseVector<DataType_>(h_size);
                data.v = new DenseVector<DataType_>(h_size);
                data.b = new DenseVector<DataType_>(h_size);
                data.temp = new DenseVector<DataType_>(h_size);

                mpi::mpi_recv(data.h->elements(), data.h->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(data.u->elements(), data.u->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(data.v->elements(), data.v->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(data.b->elements(), data.b->size(), _masterid, _masterid, _comm_cart);

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

            void _send_fringe(int target, PackedGridFringe<D2Q9> & fringe, std::list<unsigned long> & ul_buffer, std::vector<MPI_Request> & requests)
            {
                ul_buffer.push_back(fringe.h_targets->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.h_index->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_h_index->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_targets_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_index_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_h_targets->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_1->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_2->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_3->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_4->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_5->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_6->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_7->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));
                ul_buffer.push_back(fringe.external_dir_targets_8->size());
                requests.push_back(mpi::mpi_isend(&(ul_buffer.back()), 1, target, _mycartid, _comm_cart));

                requests.push_back(mpi::mpi_isend(fringe.h_targets->elements(), fringe.h_targets->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.h_index->elements(), fringe.h_index->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_h_index->elements(), fringe.external_h_index->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_1->elements(), fringe.dir_targets_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_1->elements(), fringe.dir_index_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_2->elements(), fringe.dir_targets_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_2->elements(), fringe.dir_index_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_3->elements(), fringe.dir_targets_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_3->elements(), fringe.dir_index_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_4->elements(), fringe.dir_targets_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_4->elements(), fringe.dir_index_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_5->elements(), fringe.dir_targets_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_5->elements(), fringe.dir_index_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_6->elements(), fringe.dir_targets_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_6->elements(), fringe.dir_index_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_7->elements(), fringe.dir_targets_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_7->elements(), fringe.dir_index_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_targets_8->elements(), fringe.dir_targets_8->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.dir_index_8->elements(), fringe.dir_index_8->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_1->elements(), fringe.external_dir_index_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_2->elements(), fringe.external_dir_index_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_3->elements(), fringe.external_dir_index_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_4->elements(), fringe.external_dir_index_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_5->elements(), fringe.external_dir_index_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_6->elements(), fringe.external_dir_index_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_7->elements(), fringe.external_dir_index_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_index_8->elements(), fringe.external_dir_index_8->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_h_targets->elements(), fringe.external_h_targets->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_1->elements(), fringe.external_dir_targets_1->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_2->elements(), fringe.external_dir_targets_2->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_3->elements(), fringe.external_dir_targets_3->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_4->elements(), fringe.external_dir_targets_4->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_5->elements(), fringe.external_dir_targets_5->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_6->elements(), fringe.external_dir_targets_6->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_7->elements(), fringe.external_dir_targets_7->size(), target, _mycartid, _comm_cart));
                requests.push_back(mpi::mpi_isend(fringe.external_dir_targets_8->elements(), fringe.external_dir_targets_8->size(), target, _mycartid, _comm_cart));
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

                mpi::mpi_recv(&h_targets_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&h_index_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_h_index_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_targets_8_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&dir_index_8_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_index_8_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_h_targets_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_1_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_2_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_3_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_4_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_5_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_6_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_7_size, 1, _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(&external_dir_targets_8_size, 1, _masterid, _masterid, _comm_cart);

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

                mpi::mpi_recv(fringe.h_targets->elements(), fringe.h_targets->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.h_index->elements(), fringe.h_index->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_h_index->elements(), fringe.external_h_index->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_1->elements(), fringe.dir_targets_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_1->elements(), fringe.dir_index_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_2->elements(), fringe.dir_targets_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_2->elements(), fringe.dir_index_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_3->elements(), fringe.dir_targets_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_3->elements(), fringe.dir_index_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_4->elements(), fringe.dir_targets_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_4->elements(), fringe.dir_index_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_5->elements(), fringe.dir_targets_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_5->elements(), fringe.dir_index_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_6->elements(), fringe.dir_targets_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_6->elements(), fringe.dir_index_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_7->elements(), fringe.dir_targets_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_7->elements(), fringe.dir_index_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_targets_8->elements(), fringe.dir_targets_8->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.dir_index_8->elements(), fringe.dir_index_8->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_1->elements(), fringe.external_dir_index_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_2->elements(), fringe.external_dir_index_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_3->elements(), fringe.external_dir_index_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_4->elements(), fringe.external_dir_index_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_5->elements(), fringe.external_dir_index_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_6->elements(), fringe.external_dir_index_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_7->elements(), fringe.external_dir_index_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_index_8->elements(), fringe.external_dir_index_8->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_h_targets->elements(), fringe.external_h_targets->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_1->elements(), fringe.external_dir_targets_1->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_2->elements(), fringe.external_dir_targets_2->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_3->elements(), fringe.external_dir_targets_3->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_4->elements(), fringe.external_dir_targets_4->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_5->elements(), fringe.external_dir_targets_5->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_6->elements(), fringe.external_dir_targets_6->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_7->elements(), fringe.external_dir_targets_7->size(), _masterid, _masterid, _comm_cart);
                mpi::mpi_recv(fringe.external_dir_targets_8->elements(), fringe.external_dir_targets_8->size(), _masterid, _masterid, _comm_cart);
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

                mpi::mpi_send(data.h->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_1->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_2->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_3->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_4->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_5->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_6->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_7->elements(), data.h->size(), target, _mycartid, _comm_cart);
                mpi::mpi_send(data.f_temp_8->elements(), data.h->size(), target, _mycartid, _comm_cart);

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

            void _recv_full_sync(int target, PackedGridData<D2Q9, DataType_> & data)
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

                mpi::mpi_recv(data.h->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_1->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_2->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_3->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_4->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_5->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_6->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_7->elements(), data.h->size(), target, target, _comm_cart);
                mpi::mpi_recv(data.f_temp_8->elements(), data.h->size(), target, target, _comm_cart);

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

            void _send_master_sync(int target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
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
                if (f1_size > 0) mpi::mpi_send(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, _mycartid, _comm_cart);
                unsigned long f2_offset((*fringe.dir_index_2)[0]);
                unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0)mpi::mpi_send(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, _mycartid, _comm_cart);
                unsigned long f3_offset((*fringe.dir_index_3)[0]);
                unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_send(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, _mycartid, _comm_cart);
                unsigned long f4_offset((*fringe.dir_index_4)[0]);
                unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_send(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, _mycartid, _comm_cart);
                unsigned long f5_offset((*fringe.dir_index_5)[0]);
                unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_send(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, _mycartid, _comm_cart);
                unsigned long f6_offset((*fringe.dir_index_6)[0]);
                unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_send(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, _mycartid, _comm_cart);
                unsigned long f7_offset((*fringe.dir_index_7)[0]);
                unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_send(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, _mycartid, _comm_cart);
                unsigned long f8_offset((*fringe.dir_index_8)[0]);
                unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_send(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, _mycartid, _comm_cart);

                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_send(data.h->elements() + h_offset - offset, h_size, target, _mycartid, _comm_cart);
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

            void _recv_slave_sync(int target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
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
                if (f1_size > 0) mpi::mpi_recv(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, target, _comm_cart);
                unsigned long f2_offset((*fringe.dir_index_2)[0]);
                unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0) mpi::mpi_recv(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, target, _comm_cart);
                unsigned long f3_offset((*fringe.dir_index_3)[0]);
                unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_recv(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, target, _comm_cart);
                unsigned long f4_offset((*fringe.dir_index_4)[0]);
                unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_recv(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, target, _comm_cart);
                unsigned long f5_offset((*fringe.dir_index_5)[0]);
                unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_recv(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, target, _comm_cart);
                unsigned long f6_offset((*fringe.dir_index_6)[0]);
                unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_recv(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, target, _comm_cart);
                unsigned long f7_offset((*fringe.dir_index_7)[0]);
                unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_recv(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, target, _comm_cart);
                unsigned long f8_offset((*fringe.dir_index_8)[0]);
                unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_recv(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, target, _comm_cart);

                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_recv(data.h->elements() + h_offset - offset, h_size, target, target, _comm_cart);
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

            void _send_slave_sync(int target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
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
                if (f1_size > 0) mpi::mpi_send(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, _mycartid, _comm_cart);
                unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0)mpi::mpi_send(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, _mycartid, _comm_cart);
                unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_send(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, _mycartid, _comm_cart);
                unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_send(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, _mycartid, _comm_cart);
                unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_send(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, _mycartid, _comm_cart);
                unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_send(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, _mycartid, _comm_cart);
                unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_send(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, _mycartid, _comm_cart);
                unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_send(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, _mycartid, _comm_cart);

                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_send(data.h->elements() + h_offset - offset, h_size, target, _mycartid, _comm_cart);
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

            void _recv_master_sync(int target, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, PackedGridFringe<D2Q9> & fringe)
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
                if (f1_size > 0) mpi::mpi_recv(data.f_temp_1->elements() + f1_offset - offset, f1_size, target, target, _comm_cart);
                unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                if (f2_size > 0) mpi::mpi_recv(data.f_temp_2->elements() + f2_offset - offset, f2_size, target, target, _comm_cart);
                unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                if (f3_size > 0) mpi::mpi_recv(data.f_temp_3->elements() + f3_offset - offset, f3_size, target, target, _comm_cart);
                unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                if (f4_size > 0) mpi::mpi_recv(data.f_temp_4->elements() + f4_offset - offset, f4_size, target, target, _comm_cart);
                unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                if (f5_size > 0) mpi::mpi_recv(data.f_temp_5->elements() + f5_offset - offset, f5_size, target, target, _comm_cart);
                unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                if (f6_size > 0) mpi::mpi_recv(data.f_temp_6->elements() + f6_offset - offset, f6_size, target, target, _comm_cart);
                unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                if (f7_size > 0) mpi::mpi_recv(data.f_temp_7->elements() + f7_offset - offset, f7_size, target, target, _comm_cart);
                unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                if (f8_size > 0) mpi::mpi_recv(data.f_temp_8->elements() + f8_offset - offset, f8_size, target, target, _comm_cart);

                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_size > 0) mpi::mpi_recv(data.h->elements() + h_offset - offset, h_size, target, target, _comm_cart);
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

                if (_solver_tag_value != tags::tv_gpu_cuda)
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
                }
                else
                {
#ifdef HONEI_CUDA
                    TicketVector tickets;
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
                    DetectTask detask(data.h->address(), &target);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask, _gpu_device));
                    DetectTask detask1(data.f_temp_1->address(), &target1);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask1, _gpu_device));
                    DetectTask detask2(data.f_temp_2->address(), &target2);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask2, _gpu_device));
                    DetectTask detask3(data.f_temp_3->address(), &target3);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask3, _gpu_device));
                    DetectTask detask4(data.f_temp_4->address(), &target4);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask4, _gpu_device));
                    DetectTask detask5(data.f_temp_5->address(), &target5);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask5, _gpu_device));
                    DetectTask detask6(data.f_temp_6->address(), &target6);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask6, _gpu_device));
                    DetectTask detask7(data.f_temp_7->address(), &target7);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask7, _gpu_device));
                    DetectTask detask8(data.f_temp_8->address(), &target8);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask8, _gpu_device));
                    tickets.wait();

                    unsigned long offset(info.offset);
                    unsigned long f1_offset((*fringe.dir_index_1)[0]);
                    unsigned long f1_size((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset);
                    DownloadTask dtask1((void *)((DataType_*)target1 + f1_offset - offset), (void *)((DataType_*)data.f_temp_1->address() + f1_offset - offset), f1_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask1, _gpu_device));
                    unsigned long f2_offset((*fringe.dir_index_2)[0]);
                    unsigned long f2_size((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset);
                    DownloadTask dtask2((void *)((DataType_*)target2 + f2_offset - offset), (void *)((DataType_*)data.f_temp_2->address() + f2_offset - offset), f2_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask2, _gpu_device));
                    unsigned long f3_offset((*fringe.dir_index_3)[0]);
                    unsigned long f3_size((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset);
                    DownloadTask dtask3((void *)((DataType_*)target3 + f3_offset - offset), (void *)((DataType_*)data.f_temp_3->address() + f3_offset - offset), f3_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask3, _gpu_device));
                    unsigned long f4_offset((*fringe.dir_index_4)[0]);
                    unsigned long f4_size((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset);
                    DownloadTask dtask4((void *)((DataType_*)target4 + f4_offset - offset), (void *)((DataType_*)data.f_temp_4->address() + f4_offset - offset), f4_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask4, _gpu_device));
                    unsigned long f5_offset((*fringe.dir_index_5)[0]);
                    unsigned long f5_size((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset);
                    DownloadTask dtask5((void *)((DataType_*)target5 + f5_offset - offset), (void *)((DataType_*)data.f_temp_5->address() + f5_offset - offset), f5_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask5, _gpu_device));
                    unsigned long f6_offset((*fringe.dir_index_6)[0]);
                    unsigned long f6_size((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset);
                    DownloadTask dtask6((void *)((DataType_*)target6 + f6_offset - offset), (void *)((DataType_*)data.f_temp_6->address() + f6_offset - offset), f6_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask6, _gpu_device));
                    unsigned long f7_offset((*fringe.dir_index_7)[0]);
                    unsigned long f7_size((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset);
                    DownloadTask dtask7((void *)((DataType_*)target7 + f7_offset - offset), (void *)((DataType_*)data.f_temp_7->address() + f7_offset - offset), f7_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask7, _gpu_device));
                    unsigned long f8_offset((*fringe.dir_index_8)[0]);
                    unsigned long f8_size((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset);
                    DownloadTask dtask8((void *)((DataType_*)target8 + f8_offset - offset), (void *)((DataType_*)data.f_temp_8->address() + f8_offset - offset), f8_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask8, _gpu_device));

                    for (unsigned long j(0) ; j < fringe.external_h_index->size() / 2 ; ++j)
                    {
                        unsigned long h_offset((*fringe.external_h_index)[j * 2]);
                        unsigned long h_size((*fringe.external_h_index)[j * 2 + 1] - h_offset);
                        DownloadTask dtask((void *)((DataType_*)target + h_offset - offset), (void *)((DataType_*)data.h->address() + h_offset - offset), h_size * sizeof(DataType_));
                        tickets.push_back(cuda::GPUPool::instance()->enqueue(dtask, _gpu_device));
                    }

                    tickets.wait();
#endif
                }

                std::vector<MPI_Request> requests_up;
                std::vector<MPI_Request> requests_down;
                unsigned long offset(info.offset);
                int source_up_recv, source_down_recv, target_up_send, target_down_send;
                unsigned long up_size_recv(0);
                unsigned long down_size_recv(0);
                unsigned long up_size_send(0);
                unsigned long down_size_send(0);
                MPI_Cart_shift(_comm_cart, 0, 1, &source_down_recv, &target_down_send);
                MPI_Cart_shift(_comm_cart, 0, -1, &source_up_recv, &target_up_send);


                unsigned long f1_offset_recv((*fringe.external_dir_index_1)[0]);
                unsigned long f1_size_recv((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset_recv);

                unsigned long f2_offset_recv((*fringe.external_dir_index_2)[0]);
                unsigned long f2_size_recv((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset_recv);

                unsigned long f3_offset_recv((*fringe.external_dir_index_3)[0]);
                unsigned long f3_size_recv((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset_recv);

                unsigned long f4_offset_recv((*fringe.external_dir_index_4)[0]);
                unsigned long f4_size_recv((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset_recv);

                unsigned long f5_offset_recv((*fringe.external_dir_index_5)[0]);
                unsigned long f5_size_recv((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset_recv);

                unsigned long f6_offset_recv((*fringe.external_dir_index_6)[0]);
                unsigned long f6_size_recv((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset_recv);

                unsigned long f7_offset_recv((*fringe.external_dir_index_7)[0]);
                unsigned long f7_size_recv((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset_recv);

                unsigned long f8_offset_recv((*fringe.external_dir_index_8)[0]);
                unsigned long f8_size_recv((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset_recv);

                up_size_recv += f2_size_recv + f3_size_recv + f4_size_recv + f5_size_recv;
                down_size_recv += f1_size_recv + f6_size_recv + f7_size_recv + f8_size_recv;

                unsigned long h_up_size_recv(0);
                unsigned long h_down_size_recv(0);
                unsigned long h_up_offset_recv(0);
                unsigned long h_down_offset_recv(0);
                for (unsigned long i(0) ; i < fringe.h_index->size() / 2 ; ++i)
                {
                    int h_source((*fringe.h_targets)[i]);
                    unsigned long h_offset((*fringe.h_index)[i * 2]);
                    unsigned long h_size((*fringe.h_index)[i * 2 + 1] - h_offset);
                    if (h_source < _mycartpos)
                    {
                        down_size_recv += h_size;
                        h_down_size_recv = h_size;
                        h_down_offset_recv = h_offset;
                    }
                    else
                    {
                        up_size_recv += h_size;
                        h_up_size_recv = h_size;
                        h_up_offset_recv = h_offset;
                    }
                }

                //DataType_ down_buffer_recv[down_size_recv];
                //DataType_ up_buffer_recv[up_size_recv];
                DataType_ * up_buffer_recv(TypeTraits<DataType_>::allocate(up_size_recv));
                DataType_ * down_buffer_recv(TypeTraits<DataType_>::allocate(down_size_recv));

                if (up_size_recv > 0) requests_down.push_back(mpi::mpi_irecv(up_buffer_recv, up_size_recv, source_up_recv, source_up_recv, _comm_cart));
                if (down_size_recv > 0) requests_up.push_back(mpi::mpi_irecv(down_buffer_recv, down_size_recv, source_down_recv, source_down_recv, _comm_cart));


                unsigned long f1_offset_send((*fringe.dir_index_1)[0]);
                unsigned long f1_size_send((*fringe.dir_index_1)[fringe.dir_index_1->size()-1] - f1_offset_send);

                unsigned long f2_offset_send((*fringe.dir_index_2)[0]);
                unsigned long f2_size_send((*fringe.dir_index_2)[fringe.dir_index_2->size()-1] - f2_offset_send);

                unsigned long f3_offset_send((*fringe.dir_index_3)[0]);
                unsigned long f3_size_send((*fringe.dir_index_3)[fringe.dir_index_3->size()-1] - f3_offset_send);

                unsigned long f4_offset_send((*fringe.dir_index_4)[0]);
                unsigned long f4_size_send((*fringe.dir_index_4)[fringe.dir_index_4->size()-1] - f4_offset_send);

                unsigned long f5_offset_send((*fringe.dir_index_5)[0]);
                unsigned long f5_size_send((*fringe.dir_index_5)[fringe.dir_index_5->size()-1] - f5_offset_send);

                unsigned long f6_offset_send((*fringe.dir_index_6)[0]);
                unsigned long f6_size_send((*fringe.dir_index_6)[fringe.dir_index_6->size()-1] - f6_offset_send);

                unsigned long f7_offset_send((*fringe.dir_index_7)[0]);
                unsigned long f7_size_send((*fringe.dir_index_7)[fringe.dir_index_7->size()-1] - f7_offset_send);

                unsigned long f8_offset_send((*fringe.dir_index_8)[0]);
                unsigned long f8_size_send((*fringe.dir_index_8)[fringe.dir_index_8->size()-1] - f8_offset_send);


                up_size_send += f2_size_send + f3_size_send + f4_size_send + f5_size_send;
                down_size_send += f1_size_send + f6_size_send + f7_size_send + f8_size_send;

                unsigned long h_up_size_send(0);
                unsigned long h_down_size_send(0);
                unsigned long h_up_offset_send(0);
                unsigned long h_down_offset_send(0);
                for (unsigned long i(0) ; i < fringe.external_h_index->size() / 2 ; ++i)
                {
                    int h_target((*fringe.external_h_targets)[i]);
                    unsigned long h_offset((*fringe.external_h_index)[i * 2]);
                    unsigned long h_size((*fringe.external_h_index)[i * 2 + 1] - h_offset);
                    if (h_target > _mycartpos)
                    {
                        down_size_send += h_size;
                        h_down_size_send = h_size;
                        h_down_offset_send = h_offset;
                    }
                    else
                    {
                        up_size_send += h_size;
                        h_up_size_send = h_size;
                        h_up_offset_send = h_offset;
                    }
                }

                //DataType_ up_buffer_send[up_size_send];
                //DataType_ down_buffer_send[down_size_send];
                DataType_ * up_buffer_send(TypeTraits<DataType_>::allocate(up_size_send));
                DataType_ * down_buffer_send(TypeTraits<DataType_>::allocate(down_size_send));

                if (up_size_send > 0)
                {
                    unsigned long temp_size(0);
                    TypeTraits<DataType_>::copy(data.f_temp_2->elements() + f2_offset_send - offset, up_buffer_send + temp_size, f2_size_send);
                    temp_size += f2_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_3->elements() + f3_offset_send - offset, up_buffer_send + temp_size, f3_size_send);
                    temp_size += f3_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_4->elements() + f4_offset_send - offset, up_buffer_send + temp_size, f4_size_send);
                    temp_size += f4_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_5->elements() + f5_offset_send - offset, up_buffer_send + temp_size, f5_size_send);
                    temp_size += f5_size_send;
                    TypeTraits<DataType_>::copy(data.h->elements() + h_up_offset_send - offset, up_buffer_send + temp_size, h_up_size_send);
                    requests_up.push_back(mpi::mpi_isend(up_buffer_send, up_size_send, target_up_send, _mycartid, _comm_cart));
                }

                if (down_size_send > 0)
                {
                    unsigned long temp_size(0);
                    TypeTraits<DataType_>::copy(data.f_temp_1->elements() + f1_offset_send - offset, down_buffer_send + temp_size, f1_size_send);
                    temp_size += f1_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_6->elements() + f6_offset_send - offset, down_buffer_send + temp_size, f6_size_send);
                    temp_size += f6_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_7->elements() + f7_offset_send - offset, down_buffer_send + temp_size, f7_size_send);
                    temp_size += f7_size_send;
                    TypeTraits<DataType_>::copy(data.f_temp_8->elements() + f8_offset_send - offset, down_buffer_send + temp_size, f8_size_send);
                    temp_size += f8_size_send;
                    TypeTraits<DataType_>::copy(data.h->elements() + h_down_offset_send - offset, down_buffer_send + temp_size, h_down_size_send);
                    requests_down.push_back(mpi::mpi_isend(down_buffer_send, down_size_send, target_down_send, _mycartid, _comm_cart));
                }


                TimeStamp ca, cu, cd;
                int up_fin(false);
                int down_fin(false);
                double delta_up(0);
                double delta_down(0);
                ca.take();
                while (!up_fin || !down_fin)
                {
                    if (!down_fin)
                    {
                        MPI_Testall(requests_down.size(),  &requests_down[0], &down_fin, MPI_STATUSES_IGNORE);
                        if (down_fin)
                            cd.take();
                    }
                    if (!up_fin)
                    {
                        MPI_Testall(requests_up.size(),  &requests_up[0], &up_fin, MPI_STATUSES_IGNORE);
                        if (up_fin)
                            cu.take();
                    }
                }
                //MPI_Waitall(requests_up.size(), &requests_up[0], MPI_STATUSES_IGNORE);
                //MPI_Waitall(requests_down.size(), &requests_down[0], MPI_STATUSES_IGNORE);
                delta_down = cd.total() - ca.total();
                _sync_time_down+=delta_down;
                delta_up = cu.total() - ca.total();
                _sync_time_up+=delta_up;
                _balance_load(delta_up, delta_down);

                if (up_size_recv > 0)
                {
                    unsigned long temp_size(0);
                    TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_2->elements() + f2_offset_recv - offset, f2_size_recv);
                    temp_size += f2_size_recv;
                    TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_3->elements() + f3_offset_recv - offset, f3_size_recv);
                    temp_size += f3_size_recv;
                    TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_4->elements() + f4_offset_recv - offset, f4_size_recv);
                    temp_size += f4_size_recv;
                    TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.f_temp_5->elements() + f5_offset_recv - offset, f5_size_recv);
                    temp_size += f5_size_recv;
                    TypeTraits<DataType_>::copy(up_buffer_recv + temp_size, data.h->elements() + h_up_offset_recv - offset, h_up_size_recv);
                }

                if (down_size_recv > 0)
                {
                    unsigned long temp_size(0);
                    TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_1->elements() + f1_offset_recv - offset, f1_size_recv);
                    temp_size += f1_size_recv;
                    TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_6->elements() + f6_offset_recv - offset, f6_size_recv);
                    temp_size += f6_size_recv;
                    TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_7->elements() + f7_offset_recv - offset, f7_size_recv);
                    temp_size += f7_size_recv;
                    TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.f_temp_8->elements() + f8_offset_recv - offset, f8_size_recv);
                    temp_size += f8_size_recv;
                    TypeTraits<DataType_>::copy(down_buffer_recv + temp_size, data.h->elements() + h_down_offset_recv - offset, h_down_size_recv);
                }

                TypeTraits<DataType_>::free(up_buffer_send, up_size_send);
                TypeTraits<DataType_>::free(down_buffer_send, down_size_send);
                TypeTraits<DataType_>::free(up_buffer_recv, up_size_recv);
                TypeTraits<DataType_>::free(down_buffer_recv, down_size_recv);

                if (_solver_tag_value != tags::tv_gpu_cuda)
                {
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
                else
                {
#ifdef HONEI_CUDA
                    TicketVector tickets;
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
                    DetectTask detask(data.h->address(), &target);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask, _gpu_device));
                    DetectTask detask1(data.f_temp_1->address(), &target1);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask1, _gpu_device));
                    DetectTask detask2(data.f_temp_2->address(), &target2);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask2, _gpu_device));
                    DetectTask detask3(data.f_temp_3->address(), &target3);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask3, _gpu_device));
                    DetectTask detask4(data.f_temp_4->address(), &target4);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask4, _gpu_device));
                    DetectTask detask5(data.f_temp_5->address(), &target5);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask5, _gpu_device));
                    DetectTask detask6(data.f_temp_6->address(), &target6);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask6, _gpu_device));
                    DetectTask detask7(data.f_temp_7->address(), &target7);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask7, _gpu_device));
                    DetectTask detask8(data.f_temp_8->address(), &target8);
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(detask8, _gpu_device));
                    tickets.wait();

                    unsigned long offset(info.offset);
                    unsigned long f1_offset((*fringe.external_dir_index_1)[0]);
                    unsigned long f1_size((*fringe.external_dir_index_1)[fringe.external_dir_index_1->size()-1] - f1_offset);
                    UploadTask utask1((void *)((DataType_*)data.f_temp_1->address() + f1_offset - offset), (void *)((DataType_*)target1 + f1_offset - offset), f1_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask1, _gpu_device));
                    unsigned long f2_offset((*fringe.external_dir_index_2)[0]);
                    unsigned long f2_size((*fringe.external_dir_index_2)[fringe.external_dir_index_2->size()-1] - f2_offset);
                    UploadTask utask2((void *)((DataType_*)data.f_temp_2->address() + f2_offset - offset), (void *)((DataType_*)target2 + f2_offset - offset), f2_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask2, _gpu_device));
                    unsigned long f3_offset((*fringe.external_dir_index_3)[0]);
                    unsigned long f3_size((*fringe.external_dir_index_3)[fringe.external_dir_index_3->size()-1] - f3_offset);
                    UploadTask utask3((void *)((DataType_*)data.f_temp_3->address() + f3_offset - offset), (void *)((DataType_*)target3 + f3_offset - offset), f3_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask3, _gpu_device));
                    unsigned long f4_offset((*fringe.external_dir_index_4)[0]);
                    unsigned long f4_size((*fringe.external_dir_index_4)[fringe.external_dir_index_4->size()-1] - f4_offset);
                    UploadTask utask4((void *)((DataType_*)data.f_temp_4->address() + f4_offset - offset), (void *)((DataType_*)target4 + f4_offset - offset), f4_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask4, _gpu_device));
                    unsigned long f5_offset((*fringe.external_dir_index_5)[0]);
                    unsigned long f5_size((*fringe.external_dir_index_5)[fringe.external_dir_index_5->size()-1] - f5_offset);
                    UploadTask utask5((void *)((DataType_*)data.f_temp_5->address() + f5_offset - offset), (void *)((DataType_*)target5 + f5_offset - offset), f5_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask5, _gpu_device));
                    unsigned long f6_offset((*fringe.external_dir_index_6)[0]);
                    unsigned long f6_size((*fringe.external_dir_index_6)[fringe.external_dir_index_6->size()-1] - f6_offset);
                    UploadTask utask6((void *)((DataType_*)data.f_temp_6->address() + f6_offset - offset), (void *)((DataType_*)target6 + f6_offset - offset), f6_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask6, _gpu_device));
                    unsigned long f7_offset((*fringe.external_dir_index_7)[0]);
                    unsigned long f7_size((*fringe.external_dir_index_7)[fringe.external_dir_index_7->size()-1] - f7_offset);
                    UploadTask utask7((void *)((DataType_*)data.f_temp_7->address() + f7_offset - offset), (void *)((DataType_*)target7 + f7_offset - offset), f7_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask7, _gpu_device));
                    unsigned long f8_offset((*fringe.external_dir_index_8)[0]);
                    unsigned long f8_size((*fringe.external_dir_index_8)[fringe.external_dir_index_8->size()-1] - f8_offset);
                    UploadTask utask8((void *)((DataType_*)data.f_temp_8->address() + f8_offset - offset), (void *)((DataType_*)target8 + f8_offset - offset), f8_size * sizeof(DataType_));
                    tickets.push_back(cuda::GPUPool::instance()->enqueue(utask8, _gpu_device));

                    for (unsigned long j(0) ; j < fringe.h_index->size() / 2 ; ++j)
                    {
                        unsigned long h_offset((*fringe.h_index)[j * 2]);
                        unsigned long h_size((*fringe.h_index)[j * 2 + 1] - h_offset);
                        UploadTask utask((void *)((DataType_*)data.h->address() + h_offset - offset), (void *)((DataType_*)target + h_offset - offset), h_size * sizeof(DataType_));
                        tickets.push_back(cuda::GPUPool::instance()->enqueue(utask, _gpu_device));
                    }

                    tickets.wait();
#endif
                }
            }

            void _balance_load(double delta_up, double delta_down)
            {
                std::vector<MPI_Request> requests;
                int source_up_recv, source_down_recv, target_up_send, target_down_send;

                MPI_Cart_shift(_comm_cart, 0, 1, &source_up_recv, &target_down_send);
                MPI_Cart_shift(_comm_cart, 0, -1, &source_down_recv, &target_up_send);

                // Load Balancing Requests

                // 0 means: nothing to do
                // 1 means: give me data, you are to slow
                unsigned long status_in_up(0);
                unsigned long status_in_down(0);
                if (source_up_recv != MPI_PROC_NULL) requests.push_back(mpi::mpi_irecv(&status_in_up, 1, source_up_recv, source_up_recv, _comm_cart));
                if (source_down_recv != MPI_PROC_NULL) requests.push_back(mpi::mpi_irecv(&status_in_down, 1, source_down_recv, source_down_recv, _comm_cart));

                unsigned long status_out_up(0);
                unsigned long status_out_down(0);

                if (delta_up > _sync_threshold)
                    status_out_up = 1;

                if (delta_down > _sync_threshold)
                    status_out_down = 1;

                if (target_up_send != MPI_PROC_NULL) requests.push_back(mpi::mpi_isend(&status_out_up, 1, target_up_send, _mycartid, _comm_cart));
                if (target_down_send != MPI_PROC_NULL) requests.push_back(mpi::mpi_isend(&status_out_down, 1, target_down_send, _mycartid, _comm_cart));
                MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
                requests.clear();

                // Answers

                // 0 means: go ahead, i am waiting, too | or if no request was sent, gives a dummy answer
                // 1 means: ok, i will give you some data
                // A -> B means: A -> waits for B, where A is over B

                unsigned long old_out_up = status_out_up;
                unsigned long old_out_down = status_out_down;
                // if down is waiting for us
                // W <- N
                if (status_in_down == 1 && old_out_up == 0)
                {
                    status_out_down = 1;
                    //std::cout<<_mycartid<<": down waits for me"<<std::endl;
                }
                // if up is waiting for us
                // N -> W
                if (status_in_up == 1 && old_out_down == 0)
                {
                    status_out_up = 1;
                    //std::cout<<_mycartid<<": up waits for me"<<std::endl;
                }
                // if down is waiting for us but we are waiting for up, too
                // N <- W <- N
                if (status_in_down == 1 && old_out_up == 1)
                {
                    status_out_down = 0;
                    //std::cout<<_mycartid<<": down waits for me, but i wait for up"<<std::endl;
                }
                // if up is waiting for us but we are waiting for down, too
                // N -> W -> N
                if (status_in_up == 1 && old_out_down == 1)
                {
                    status_out_up = 0;
                    //std::cout<<_mycartid<<": up waits for me, but i wait for down"<<std::endl;
                }
                // if down is waiting for us and we are waiting for up, but up is waiting for us too
                // N <-> W <- N
                if (status_in_down == 1 && old_out_up == 1 && status_in_up == 1)
                {
                    status_out_down = 1;
                    //std::cout<<_mycartid<<": down and up wait for me, i wait for up"<<std::endl;
                }
                // if up is waiting for us and we are waiting for down, but down is waiting for us too
                // N -> W <-> N
                if (status_in_up == 1 && old_out_down == 1 && status_in_down == 1)
                {
                    status_out_up = 1;
                    //std::cout<<_mycartid<<": up and down wait for me, i wait for down"<<std::endl;
                }
                /* // if down and up are waiting for us
                if (status_in_up == 1 && status_in_down == 1)
                {
                    status_out_up = 0;
                    status_out_down = 1;
                    std::cout<<_mycartid<<": up and down waits for me, favouring down"<<std::endl;
                }*/
                // down wants nothing from us or we synched recently
                if (status_in_down == 0 || _recently_synched)
                {
                    status_out_down = 0;
                    status_in_down = 0;
                }
                // up wants nothing from us or we synched recently
                if (status_in_up == 0 || _recently_synched)
                {
                    status_out_up = 0;
                    status_in_up = 0;
                }
                // if down waits for us and we wait for down
                // W <-> N
                if (status_in_down == 1 && old_out_down == 1 && !_recently_synched)
                {
                    status_out_down = 0;
                    //status_out_up = 1; //experimentell
                    //std::cout<<_mycartid<<": down waits for me, i wait for down"<<std::endl;
                }
                // if up waits for us and we wait for up
                // N <-> W
                if (status_in_up == 1 && old_out_up == 1 && !_recently_synched)
                {
                    status_out_up = 0;
                    //status_out_down = 1; //experimentell
                    //std::cout<<_mycartid<<": up waits for me, i wait for up"<<std::endl;
                }
                _recently_synched = false;

                if (source_up_recv != MPI_PROC_NULL) requests.push_back(mpi::mpi_irecv(&status_in_up, 1, source_up_recv, source_up_recv, _comm_cart));
                if (source_down_recv != MPI_PROC_NULL) requests.push_back(mpi::mpi_irecv(&status_in_down, 1, source_down_recv, source_down_recv, _comm_cart));
                if (target_up_send != MPI_PROC_NULL) requests.push_back(mpi::mpi_isend(&status_out_up, 1, target_up_send, _mycartid, _comm_cart));
                if (target_down_send != MPI_PROC_NULL) requests.push_back(mpi::mpi_isend(&status_out_down, 1, target_down_send, _mycartid, _comm_cart));

                MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
                requests.clear();
                //std::cout<<_mycartid<<": "<<delta_down<<", delta up: "<<delta_up<<std::endl;

                // status_in == 0 : you get no data
                // status_in == 1 : you get data
                // if status_out up/down == 1 -> send/recv data
                // \TODO recv/send real data

                // receive data from up
                if (status_in_up == 1)
                {
                    _recently_synched = true;
                    //std::cout<<_mycartid<<": up gives me data"<<std::endl;
                }
                // receive data from down
                if (status_in_down == 1)
                {
                    _recently_synched = true;
                    //std::cout<<_mycartid<<": down gives me data"<<std::endl;
                }
                // send data to up
                if (status_out_up == 1)
                {
                    _recently_synched = true;
                    //std::cout<<_mycartid<<": i give up data"<<std::endl;
                }
                // send data to down
                if (status_out_down == 1)
                {
                    _recently_synched = true;
                    //std::cout<<_mycartid<<": i give down data"<<std::endl;
                }
            }

            void _read_config(std::string filename, unsigned long & scenario, std::vector<std::string> & backends,
                    std::vector<double> & fractions, bool & file_output, std::string & base_filename)
            {
                if (_mycartid == _masterid)
                    std::cout<<"reading..."<<filename<<std::endl;

                std::ifstream file(filename.c_str());
                if (!file.is_open())
                    throw honei::InternalError("Unable to open mpi config file: " + filename);

                while(!file.eof())
                {
                    std::string line;
                    std::getline(file, line);
                    if (file.eof())
                        break;

                    //skip blank lines
                    if (line.size() == 0)
                        continue;

                    //skip comment lines
                    if(line.find("#") == 0)
                        continue;

                    if (_mycartid == _masterid)
                        std::cout<<line<<std::endl;

                    //read in valueable lines
                    std::vector<std::string> line_parts(string_tokenizer(line, " "));

                    if (line_parts.at(0).compare("scenario") == 0)
                    {
                        if (line_parts.size() != 2)
                            throw InternalError("Wrong argument count to scenario statement: " + line);

                        scenario = atoi(line_parts.at(1).c_str());
                    }
                    else if (line_parts.at(0).compare("proc") == 0)
                    {
                        if (line_parts.size() != 3)
                            throw InternalError("Wrong argument count to proc statement: " + line);

                        backends.push_back(line_parts.at(1));
                        fractions.push_back(atof(line_parts.at(2).c_str()));
                    }
                    else if (line_parts.at(0).compare("fileoutput") == 0)
                    {
                        if (line_parts.size() != 2)
                            throw InternalError("Wrong argument count to fileoutput statement: " + line);

                        base_filename = line_parts.at(1);
                        file_output = true;
                    }
                    else if (line_parts.at(0).compare("syncthreshold") == 0)
                    {
                        if (line_parts.size() != 2)
                            throw InternalError("Wrong argument count to syncthreshold statement: " + line);

                        _sync_threshold = atof(line_parts.at(1).c_str());
                    }
                    else
                        throw InternalError("Error: config file entry not known: " + line);
                }
                file.close();
            }

            void _init_solver(SolverLBMGridBase *& solver, unsigned long id, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DataType_> & data, DataType_ d_x, DataType_ d_y, DataType_ d_t, DataType_ tau)
            {
                if(_backends.at(id).compare(tags::CPU::SSE::name) == 0)
                {
#ifdef HONEI_SSE
                    _device_name = tags::CPU::SSE::name;
                    _solver_tag_value = tags::CPU::SSE::tag_value;
                    solver = new SolverLBMGrid<tags::CPU::SSE, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (&info, &data, d_x, d_y, d_t, tau);
#else
                    throw InternalError("Backend not activated: " + _backends.at(id));
#endif
                }
                else if(_backends.at(id).compare(tags::GPU::CUDA::name) == 0)
                {
#ifdef HONEI_CUDA
                    _device_name = tags::GPU::CUDA::name;
                    _solver_tag_value = tags::GPU::CUDA::tag_value;
                    // die wievielte gpu sind wir auf dem knoten?
                    unsigned long gpu_number(0);
                    for (unsigned long i(0) ; i < id ; ++i)
                        if (_backends.at(i).compare(tags::GPU::CUDA::name) == 0)
                            gpu_number++;
                    cuda::GPUPool::instance()->single_start(gpu_number);
                    solver = new SolverLBMGrid<tags::GPU::CUDA, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> (&info, &data, d_x, d_y, d_t, tau);
#else
                    throw InternalError("Backend not activated: " + _backends.at(id));
#endif
                }
                else
                {
                    throw InternalError("Backend not known: " + _backends.at(id));
                }
            }
    };
}
#endif
