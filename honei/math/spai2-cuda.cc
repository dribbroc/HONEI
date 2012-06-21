/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/spai2.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/backends/cuda/transfer.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>
#include <honei/util/profiler.hh>

#include <iostream>
#include <honei/util/time_stamp.hh>

using namespace honei;

namespace
{
    class cudaSpai2double
    {
        private:
            SparseMatrix<double> & M;
            const SparseMatrix<double> & A;
            unsigned long blocksize;
        public:
            cudaSpai2double(SparseMatrix<double> & M, const SparseMatrix<double> & A, unsigned long blocksize) :
                M(M),
                A(A),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                TimeStamp at, bt;
                at.take();
                M._synch_column_vectors();

                DenseVector<double> a_elements(A.used_elements());
                double * a_elements_e(a_elements.elements());
                DenseVector<unsigned long> a_indices(A.used_elements());
                unsigned long * a_indices_e(a_indices.elements());
                DenseVector<double> m_elements(M.used_elements());
                double * m_elements_e(m_elements.elements());
                DenseVector<unsigned long> columns(A.columns() + 1);
                unsigned long offset(0);


                columns[0] = 0;
                for (unsigned long column(0) ; column < A.columns() ; ++column)
                {
                    double * elements(A.column(column).elements());
                    unsigned long * indices(A.column(column).indices());

                    const unsigned long size(A.column(column).used_elements());
                    for (unsigned long i(0) ; i < size ; ++i)
                    {
                        a_elements_e[i + offset] = elements[i];
                        a_indices_e[i + offset] = indices[i];
                    }
                    offset += size;
                    columns[column + 1] = offset;
                }

                void * a_elements_gpu(a_elements.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * a_indices_gpu(a_indices.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * m_elements_gpu(m_elements.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * columns_gpu(columns.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_spai2_double(columns_gpu, m_elements_gpu, a_elements_gpu, a_indices_gpu, A.columns(), blocksize);

                a_elements.unlock(lm_read_only);
                a_indices.unlock(lm_read_only);
                m_elements.unlock(lm_write_only);
                columns.unlock(lm_read_only);

                m_elements.lock(lm_read_only);
                offset = 0;
                for (unsigned long column(0) ; column < M.columns() ; ++column)
                {
                    double * elements(M.column(column).elements());

                    const unsigned long size(M.column(column).used_elements());
                    for (unsigned long i(0) ; i < size ; ++i)
                    {
                        elements[i] = m_elements_e[i + offset];
                    }
                    offset += size;
                }

                M._synch_row_vectors();
                m_elements.unlock(lm_read_only);
                bt.take();
                std::cout<<"TOE GPU: "<<bt.total()-at.total()<<std::endl;
                std::cout<<"error code: "<<m_elements_e[0]<<std::endl;
            }
    };
}

#ifdef HONEI_CUDA_DOUBLE
SparseMatrix<double> & SPAI2<tags::GPU::CUDA>::value(SparseMatrix<double> & M, const SparseMatrix<double> & A)
{
    CONTEXT("When calculating SPAI2<double> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::spai2_double", 128ul));

    M = A.copy();

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaSpai2double task(M, A, blocksize);
        task();
    }
    else
    {
        cudaSpai2double task(M, A, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }



    return M;
}
#endif
