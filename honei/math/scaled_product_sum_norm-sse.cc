/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. MATH is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * MATH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/sse/operations.hh>
#include <honei/math/scaled_product_sum_norm.hh>
#include <honei/la/product.hh>

namespace honei
{
    float ScaledProductSumNorm_TUTORIAL<tags::CPU::SSE>::value(float a, DenseVector<float> & y, float b, BandedMatrixQ1<float> & A, DenseVector<float> & x)
    {
        x.lock(lm_read_only);
        y.lock(lm_read_only);
        A.lock(lm_read_only);

        //Still use HONEIs BandedMatrix-DenseVector product:
        DenseVector<float> A_x(Product<tags::CPU::SSE>::value(A, x));

        //do not care about alignment, HONEI containers provide aligned data
        float * A_x_data = A_x.elements();
        float * y_data = y.elements();

        //redirect the relevant data to the SSE backend
        float result(honei::sse::scaled_product_sum_norm(x.size(), a, y_data, b, A_x_data));

        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
        A.unlock(lm_read_only);

        return result;

    }
    double ScaledProductSumNorm_TUTORIAL<tags::CPU::SSE>::value(double a, DenseVector<double> & y, double b, BandedMatrixQ1<double> & A, DenseVector<double> & x)
    {

        x.lock(lm_read_only);
        y.lock(lm_read_only);
        A.lock(lm_read_only);

        //Still use HONEIs BandedMatrix-DenseVector product:
        DenseVector<double> A_x(Product<tags::CPU::SSE>::value(A, x));

        //do not care about alignment, HONEI containers provide aligned data
        double * A_x_data = A_x.elements();
        double * y_data = y.elements();

        //redirect the relevant data to the SSE backend
        double result(honei::sse::scaled_product_sum_norm(x.size(), a, y_data, b, A_x_data));

        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
        A.unlock(lm_read_only);

        return result;
    }
}
