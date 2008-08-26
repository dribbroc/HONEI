/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/prolongation.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

DenseVector<float> & Prolongation<tags::GPU::CUDA>::value(DenseVector<float> & fine,
        const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask)
{
    CONTEXT("When prolongating from coarse to fine (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::prolongation_float", 64ul));

    void * fine_gpu (fine.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * coarse_gpu (coarse.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_prolongation_float(fine_gpu, fine.size(), coarse_gpu, coarse.size(), mask.elements(), blocksize);
    coarse.unlock(lm_read_only);
    fine.unlock(lm_read_and_write);

    return fine;
}
