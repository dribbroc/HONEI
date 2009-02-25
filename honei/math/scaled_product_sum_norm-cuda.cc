/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c)  2009 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/scaled_product_sum_norm.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

float ScaledProductSumNorm_TUTORIAL<tags::GPU::CUDA>::value(float a, DenseVector<float> & y, float b, BandedMatrixQ1<float> & A, DenseVector<float> & x)
{
        //Still use HONEIs BandedMatrix-DenseVector product:
        DenseVector<float> A_x(Product<tags::GPU::CUDA>::value(A, x));

        unsigned long blocksize(Configuration::instance()->get_value("cuda::spsn_float", 128ul));
        unsigned long gridsize(Configuration::instance()->get_value("cuda::spsn_float_grid", 16ul));
        void * A_x_gpu (A_x.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));

        //redirect the relevant data to the CUDA backend
        float result(cuda_scaled_product_sum_norm_float_tut(x.size(), a, y_gpu, b, A_x_gpu, blocksize, gridsize));

        y.unlock(lm_read_only);
        A_x.unlock(lm_read_only);

        return result;
}

