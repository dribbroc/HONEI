/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

#include <string>
#endif

#include <honei/lbm/grid_partitioner.hh>
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/swe/volume.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class GridPartitionerSynchBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
        int _parts;
    public:
        GridPartitionerSynchBench(const std::string & id, unsigned long size, int count, int parts) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
            _parts = parts;
        }

        virtual void run()
        {
            unsigned long g_h(_size);
            unsigned long g_w(_size);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.04), 15, 15);
            b1.value();

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            std::vector<PackedGridInfo<D2Q9> > info_list;
            std::vector<PackedGridData<D2Q9, DataType_> > data_list;
            std::vector<PackedGridFringe<D2Q9> > fringe_list;
            GridPartitioner<D2Q9, DataType_>::decompose(_parts, info, data, info_list, data_list, fringe_list);
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK((GridPartitioner<D2Q9, DataType_>::synch(info, data, info_list, data_list, fringe_list)));
            }
            evaluate();
        }
};

GridPartitionerSynchBench<tags::CPU, float> partitioner_synch_bench_float_1("LBM Grid partitioner synch Benchmark - size: 500x500, float, 250 parts", 500, 25, 250);
GridPartitionerSynchBench<tags::CPU, double> partitioner_synch_bench_double_1("LBM Grid partitioner synch Benchmark - size: 500x500, double, 250 parts", 500, 25, 250);

template <typename Tag_, typename DataType_>
class GridPartitionerDecomposeBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
        int _parts;
    public:
        GridPartitionerDecomposeBench(const std::string & id, unsigned long size, int count, int parts) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _size = size;
        _count = count;
        _parts = parts;
    }

        virtual void run()
        {

            for(int i = 0; i < _count; ++i)
            {
            unsigned long g_h(_size);
            unsigned long g_w(_size);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.04), 15, 15);
            b1.value();

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            std::vector<PackedGridInfo<D2Q9> > info_list;
            std::vector<PackedGridData<D2Q9, DataType_> > data_list;
            std::vector<PackedGridFringe<D2Q9> > fringe_list;
                BENCHMARK((GridPartitioner<D2Q9, DataType_>::decompose(_parts, info, data, info_list, data_list, fringe_list)));
            }
            evaluate();
        }
};

GridPartitionerDecomposeBench<tags::CPU, float> partitioner_decompose_bench_float_1("LBM Grid partitioner decompose Benchmark - size: 500x500, float, 25 parts", 500, 5, 25);
GridPartitionerDecomposeBench<tags::CPU, double> partitioner_decompose_bench_double_1("LBM Grid partitioner decompose Benchmark - size: 500x500, double, 25 parts", 500, 5, 25);
