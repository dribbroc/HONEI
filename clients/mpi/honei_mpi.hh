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

#include <mpi.h>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/lbm/solver_labswe_grid.hh>
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

#include <iostream>

namespace honei
{
    template <typename Tag_, typename DataType_>
    class MPISolver
    {
        private:
            int _numprocs;
            int _myid;

        public:
            MPISolver(int argc, char **argv)
            {
                mpi::mpi_init(&argc, &argv);
                mpi::mpi_comm_size(&_numprocs);
                mpi::mpi_comm_rank(&_myid);

                if (_myid == 0)
                {
                    _master();
                }
                else
                {
                    _slave();
                }
            }

            ~MPISolver()
            {
                mpi::mpi_finalize();
            }

        private:
            void _master()
            {
                unsigned long g_h(50);
                unsigned long g_w(50);
                unsigned long timesteps(100);


                DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
                Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
                c1.value();

                DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
                DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
                DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
                DenseMatrix<DataType_> b_x(PartialDerivative<Tag_, X, CENTRALDIFF>::value(b , DataType_(1)));
                DenseMatrix<DataType_> b_y(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(b , DataType_(1)));

                Grid<D2Q9, DataType_> grid;
                DenseMatrix<bool> obstacles(g_h, g_w, false);
                grid.obstacles = &obstacles;
                grid.h = &h;
                grid.u = &u;
                grid.v = &v;
                grid.b_x = &b_x;
                grid.b_y = &b_y;
                PackedGridData<D2Q9, DataType_>  data;
                PackedGridInfo<D2Q9> info;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
                std::vector<PackedGridInfo<D2Q9> > info_list;
                std::vector<PackedGridData<D2Q9, DataType_> > data_list;
                std::vector<PackedGridFringe<D2Q9> > fringe_list;
                GridPartitioner<D2Q9, DataType_>::decompose(_numprocs - 1, info, data, info_list, data_list, fringe_list);

                mpi::mpi_bcast(&timesteps, 1, 0);
                for (unsigned long target(1) ; target < _numprocs ; ++target)
                {
                    _send_info(target, info_list[target - 1]);
                    _send_data(target, data_list[target - 1]);
                }


                for (unsigned long target(1) ; target < _numprocs ; ++target)
                {
                    _recv_sync(target, data_list[target - 1]);
                }

                GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list);

                for (unsigned long target(1) ; target < _numprocs ; ++target)
                {
                    _send_sync(target, data_list[target - 1]);
                }

                // Preproc finished
                // start timesteps
                for(unsigned long i(0); i < timesteps; ++i)
                {
                    TimeStamp at, bt;
                    at.take();
                    //here are the solvers solving...
                    //and finished
                    for (unsigned long target(1) ; target < _numprocs ; ++target)
                    {
                        _recv_sync(target, data_list[target - 1]);
                    }

                    GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list);
                    GridPartitioner<D2Q9, DataType_>::compose(info, data, info_list, data_list);
                    GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                    PostProcessing<output_types::GNUPLOT>::value(h, 1, g_w, g_h, i);

                    for (unsigned long target(1) ; target < _numprocs ; ++target)
                    {
                        _send_sync(target, data_list[target - 1]);
                    }
                    bt.take();
                    std::cout<<"Timestep: " << i << "/" << timesteps << " TOE: "<<bt.sec() - at.sec()<<" "<<bt.usec() - at.usec()<<std::endl;
                }
            }

            void _slave()
            {
                PackedGridData<D2Q9, DataType_> data;
                PackedGridInfo<D2Q9> info;
                unsigned long timesteps;

                mpi::mpi_bcast(&timesteps, 1, 0);

                _recv_info(info);
                _recv_data(data);

                SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1.);

                solver.do_preprocessing();

                _send_sync(0, data);
                _recv_sync(0, data);

                for(unsigned long i(0); i < timesteps; ++i)
                {
                    solver.solve();

                    _send_sync(0, data);
                    _recv_sync(0, data);
                }
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
                mpi::mpi_send(data.b_x->elements(), data.b_x->size(), target, _myid);
                mpi::mpi_send(data.b_y->elements(), data.b_y->size(), target, _myid);
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
                data.b_x = new DenseVector<DataType_>(h_size);
                data.b_y = new DenseVector<DataType_>(h_size);

                mpi::mpi_recv(data.h->elements(), data.h->size(), 0, 0);
                mpi::mpi_recv(data.u->elements(), data.u->size(), 0, 0);
                mpi::mpi_recv(data.v->elements(), data.v->size(), 0, 0);
                mpi::mpi_recv(data.b_x->elements(), data.b_x->size(), 0, 0);
                mpi::mpi_recv(data.b_y->elements(), data.b_y->size(), 0, 0);

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

            void _send_sync(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                mpi::mpi_send(data.f_temp_1->elements(), data.f_temp_1->size(), target, _myid);
                mpi::mpi_send(data.f_temp_2->elements(), data.f_temp_2->size(), target, _myid);
                mpi::mpi_send(data.f_temp_3->elements(), data.f_temp_3->size(), target, _myid);
                mpi::mpi_send(data.f_temp_4->elements(), data.f_temp_4->size(), target, _myid);
                mpi::mpi_send(data.f_temp_5->elements(), data.f_temp_5->size(), target, _myid);
                mpi::mpi_send(data.f_temp_6->elements(), data.f_temp_6->size(), target, _myid);
                mpi::mpi_send(data.f_temp_7->elements(), data.f_temp_7->size(), target, _myid);
                mpi::mpi_send(data.f_temp_8->elements(), data.f_temp_8->size(), target, _myid);

                mpi::mpi_send(data.h->elements(), data.h->size(), target, _myid);
            }

            void _recv_sync(unsigned long target, PackedGridData<D2Q9, DataType_> & data)
            {
                mpi::mpi_recv(data.f_temp_1->elements(), data.f_temp_1->size(), target, target);
                mpi::mpi_recv(data.f_temp_2->elements(), data.f_temp_2->size(), target, target);
                mpi::mpi_recv(data.f_temp_3->elements(), data.f_temp_3->size(), target, target);
                mpi::mpi_recv(data.f_temp_4->elements(), data.f_temp_4->size(), target, target);
                mpi::mpi_recv(data.f_temp_5->elements(), data.f_temp_5->size(), target, target);
                mpi::mpi_recv(data.f_temp_6->elements(), data.f_temp_6->size(), target, target);
                mpi::mpi_recv(data.f_temp_7->elements(), data.f_temp_7->size(), target, target);
                mpi::mpi_recv(data.f_temp_8->elements(), data.f_temp_8->size(), target, target);

                mpi::mpi_recv(data.h->elements(), data.h->size(), target, target);
            }
    };
}
#endif
