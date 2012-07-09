/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/backends/cuda/cuda_util.hh>


// TODO sparse eingabe arrays und heap coalesced zugreifen, nicht alles auf einem block!

namespace honei
{
    namespace cuda
    {
        template <typename DT_>
        __device__ DT_ * get_next_heap(char * heap, unsigned long & offset, unsigned long count)
        {
            unsigned long old_offset = offset;
            offset += count * sizeof(DT_) * gridDim.x * blockDim.x;
            return (DT_ *)(heap + old_offset);
        }

        template <typename DT_>
        class Matrix
        {
            public:
                DT_ * data;
                unsigned long rows;
                unsigned long columns;

                __device__ Matrix(DT_ * heap, unsigned long rows, unsigned long columns) :
                    data(heap),
                    rows(rows),
                    columns(columns)
                {
                }

                inline __device__ DT_ & operator() (unsigned long row, unsigned long column)
                {
                    unsigned long threads = gridDim.x * blockDim.x;
                    return data[(columns * row + column) * threads];
                }
        };

        template <typename DT_>
        class Vector
        {
            public:
                DT_ * data;
                unsigned long size;

                __device__ Vector(DT_ * heap, unsigned long size) :
                    data(heap),
                    size(size)
                {
                }

                inline __device__ DT_ & operator[] (unsigned long index)
                {
                    unsigned long threads = gridDim.x * blockDim.x;
                    return data[index * threads];
                }
        };

        template <typename DT_>
        __global__ void spai2_gpu(unsigned long * column_ptr, DT_ * m_elements,
                DT_ * a_elements, unsigned long * a_indices,
                unsigned long columns, char * heap_gpu, unsigned long heap_size, unsigned long idx_offset)
        {
            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;
            unsigned long real_idx = idx;
            idx += idx_offset;

            if (idx >= columns)
                return;

            unsigned long heap_offset(0);
            char * local_heap = heap_gpu + ((real_idx) * sizeof(unsigned long));

            // ASSEMBLY START
            const unsigned long n2(column_ptr[idx + 1] - column_ptr[idx]);
            if (n2 == 0)
                return;

            Vector<unsigned long> J(get_next_heap<unsigned long>(local_heap, heap_offset, n2), n2);

            for (unsigned long i(0) ; i < n2 ; ++i)
            {
                J[i] = a_indices[column_ptr[idx] + i];
            }

            unsigned long n1(0);

            for (unsigned long i(0) ; i < n2 ; ++i)
            {
                n1 += column_ptr[J[i] + 1] - column_ptr[J[i]];
            }

            if (n1 == 0)
                return;

            Vector<DT_> pro_v(get_next_heap<DT_>(local_heap, heap_offset, n2), n2);
            Matrix<DT_> product(get_next_heap<DT_>(local_heap, heap_offset, n2 * n2), n2, n2);
            unsigned long offset_behind_product(heap_offset);
            Vector<unsigned long> I(get_next_heap<unsigned long>(local_heap, heap_offset, n1), n1);
            Vector<DT_> et(get_next_heap<DT_>(local_heap, heap_offset, n1), n1);
            Matrix<DT_> At(get_next_heap<DT_>(local_heap, heap_offset, n1 * n2), n1, n2);
            Matrix<DT_> Atrans(get_next_heap<DT_>(local_heap, heap_offset, n2 * n1), n2, n1);

            for (unsigned long i(0) ; i < n1 ; ++i)
            {
                et[i] = DT_(0);
            }

            unsigned long tmp(0);
            for (unsigned long i(0) ; i < n2 ; ++i)
            {
                unsigned long A_col_ji_ue(column_ptr[J[i] + 1] - column_ptr[J[i]]);
                unsigned long * A_col_ji_index(a_indices + column_ptr[J[i]]);
                for (unsigned long j(0) ; j < A_col_ji_ue ; ++j)
                {
                    I[tmp] = A_col_ji_index[j];
                    et[tmp] = (I[tmp] == idx);
                    ++tmp;
                }
            }

            for (unsigned long j(0) ; j < n2 ; ++j)
            {
                unsigned long * indices(a_indices + column_ptr[J[j]]);
                DT_ * elements(a_elements + column_ptr[J[j]]);
                unsigned long used_elements(column_ptr[J[j] + 1] - column_ptr[J[j]]);
                for (unsigned long i(0) ; i < n1 ; ++i)
                {
                    unsigned long index(I[i]);
                    unsigned long it(0);
                    while (it < used_elements)
                    {
                        if (indices[it] >= index)
                            break;
                        ++it;
                    }
                    if (it < used_elements)
                    {
                        At(i, j) = (indices[it] == index) ? elements[it] : DT_(0);
                    }
                    else
                        At(i, j) = DT_(0);
                }
            }
            // ASSEMBLY END

            for (unsigned long i(0) ; i < At.rows ; ++i)
            {
                for (unsigned long j(0) ; j < At.columns ; ++j)
                {
                    Atrans(j, i) = At(i, j);
                }
            }

            // TODO matrix produkt cache blocken?
            for (unsigned long i(0) ; i < product.rows ; ++i)
            {
                for (unsigned long j(0) ; j < product.columns ; ++j)
                {
                    DT_ temp(0);
                    for (unsigned long k(0) ; k < Atrans.columns ; ++k)
                    {
                        temp += Atrans(i, k) * At(k, j);
                    }
                    product(i, j) = temp;
                }
            }

            for (unsigned long i(0) ; i < Atrans.rows ; ++i)
            {
                DT_ temp(0);
                for (unsigned long j(0) ; j < Atrans.columns ; ++j)
                {
                    temp += Atrans(i, j) * et[j];
                }
                pro_v[i] = temp;
            }

            heap_offset = offset_behind_product;

            // LU DECOMPOSITION START
            Vector<DT_> res(get_next_heap<DT_>(local_heap, heap_offset, product.columns), product.columns);
            Matrix<DT_> u(product);
            for (unsigned long i(0) ; i < product.rows ; ++i)
            {
                for (unsigned long j(0) ; j < product.columns ; ++j)
                {
                    u(i, j) = product(i, j);
                }
            }
            Matrix<DT_> l(get_next_heap<DT_>(local_heap, heap_offset, product.rows * product.columns), product.rows, product.columns);
            for (unsigned long i(0) ; i < l.rows ; ++i)
            {
                for (unsigned long j(0) ; j < l.columns ; ++j)
                {
                    l(i, j) = DT_(0);
                }
            }
            for (unsigned long i(0) ; i < l.rows ; ++i)
            {
                l(i, i) = DT_(1);
            }

            for (unsigned long k(0) ; k < u.rows - 1 ; ++k)
            {
                for (unsigned long j(k + 1) ; j < u.rows ; ++j)
                {
                    l(j, k) = u(j, k) / u(k,k);
                    for (unsigned long i(k) ; i < u.rows ; ++i)
                    {
                        u(j, i) = u(j, i) - l(j, k) * u(k, i);
                    }
                }
            }

            for (unsigned long i(0) ; i < product.columns ; ++i)
            {
                DT_ sum(0);
                for (unsigned long j(0) ; j < i ; ++j)
                {
                    sum += l(i, j) * res[j];
                }
                res[i] = DT_(1) / l(i, i) * (pro_v[i] - sum);
            }

            for (long i(product.columns - 1) ; i >= 0 ; --i)
            {
                DT_ sum(0);
                for (unsigned long j(i+1) ; j < product.columns ; ++j)
                {
                    sum += u(i, j) * res[j];
                }
                res[i] = DT_(1) / u(i, i) * (res[i] - sum);
            }
            // LU DECOMPOSITION END

            for (unsigned long i(0) ; i < n2 ; ++i)
            {
                m_elements[column_ptr[idx] + i] = res[i];
            }
        }
    }

    template <typename DT_>
        void cuda_spai2(void * column_ptr, void * m_elements,
                void * a_elements, void * a_indices, unsigned long columns, unsigned long blocksize)
    {
        unsigned long heap_size(1800000000ul);
        char * heap_gpu(NULL);
        cudaMalloc((void**)&heap_gpu, heap_size);

        unsigned long * column_ptr_gpu((unsigned long *)column_ptr);
        DT_ * m_elements_gpu((DT_ *)m_elements);
        DT_ * a_elements_gpu((DT_ *)a_elements);
        unsigned long * a_indices_gpu((unsigned long *)a_indices);

        unsigned long parts(1);
        unsigned long column_offset(0);
        for (unsigned long i(0) ; i < parts ; ++i)
        {
            unsigned long part_size(columns / parts);
            if (i == 0)
                part_size = part_size + columns % parts;

            dim3 grid;
            dim3 block;
            block.x = blocksize;
            grid.x = (unsigned)ceil(part_size/(double)(block.x));

            honei::cuda::spai2_gpu<<<grid, block>>>(column_ptr_gpu, m_elements_gpu, a_elements_gpu, a_indices_gpu,
                part_size + column_offset, heap_gpu, heap_size, column_offset);

            column_offset += part_size;
        }

        cudaFree(heap_gpu);

        CUDA_ERROR();
    }
}

extern "C" void cuda_spai2_float(void * column_ptr, void * m_elements,
        void * a_elements, void * a_indices, unsigned long columns, unsigned long blocksize)
{
    honei::cuda_spai2<float>(column_ptr, m_elements, a_elements, a_indices, columns, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_spai2_double(void * column_ptr, void * m_elements,
        void * a_elements, void * a_indices, unsigned long columns, unsigned long blocksize)
{
    honei::cuda_spai2<double>(column_ptr, m_elements, a_elements, a_indices, columns, blocksize);
}
#endif
