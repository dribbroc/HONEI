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

#include <benchmark/benchmark.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/shared_array.hh>
#include <honei/la/dense_vector.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class MemoryArbiterBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
        int _arrays;

    public:
        MemoryArbiterBench(const std::string & id, unsigned long size, int count, int arrays) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
            _arrays = arrays;
        }



        virtual void run()
        {
            DenseVector<DataType_> * vectors[_arrays];
            for (unsigned long i(0) ; i < _arrays ; ++i)
            {
                vectors[i] = new DenseVector<DataType_>(1);
                MemoryArbiter::instance()->read<Tag_>(vectors[i]->memid(), vectors[i]->address(), vectors[i]->size());
            }
            DenseVector<DataType_> dv(_size, DataType_(42));
            for(int i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10000 ; ++j)
                        {
                            MemoryArbiter::instance()->write<Tag_>(dv.memid(), dv.address(), dv.size());
                            MemoryArbiter::instance()->release_write<Tag_>(dv.memid());
                            MemoryArbiter::instance()->read<Tag_>(dv.memid(), dv.address(), dv.size());
                            MemoryArbiter::instance()->release_read<Tag_>(dv.memid());
                        }
                        );
            }
            evaluate();
            for (unsigned long i(0) ; i < _arrays ; ++i)
            {
                MemoryArbiter::instance()->release_read<Tag_>(vectors[i]->memid());
                delete vectors[i];
            }
    }
};
MemoryArbiterBench<tags::CPU, float>  MABenchfloat ("CPU Memory Arbiter Benchmark 1 array - vector size: 64^4, float",  64ul*64*64*64, 50, 1);
MemoryArbiterBench<tags::GPU::CUDA, float>  cudaMABenchfloat ("CUDA Memory Arbiter Benchmark 1 array - vector size: 64^4, float",  64ul*64*64*64, 50, 1);
MemoryArbiterBench<tags::CPU, float>  MABenchfloat2 ("CPU Memory Arbiter Benchmark 1000 arrays - vector size: 64^4, float",  64ul*64*64*64, 50, 1000);
MemoryArbiterBench<tags::GPU::CUDA, float>  cudaMABenchfloat2 ("CUDA Memory Arbiter Benchmark 1000 arrays - vector size: 64^4, float",  64ul*64*64*64, 50, 1000);
