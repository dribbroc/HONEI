
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

#ifndef CUDA_GUARD_OPERATIONS_HH
#define CUDA_GUARD_OPERATIONS_HH 1

extern "C"
{
        /*void difference(float * a, const float * b, unsigned long size);
        void difference(double * a, const double * b, unsigned long size);

        float dot_product(const float * a, float * b, unsigned long size);
        double dot_product(double * a, double * b, unsigned long size);
*/
        void cuda_element_inverse_one_float(float * x, unsigned long size);

        void cuda_element_product_two_float(float * a, const float * b, unsigned long size);
/*
        float norm_l2(const float * a, unsigned long size);
        double norm_l2(double * a, unsigned long size);

        float reduction_sum(const float * a, unsigned long size);
        double reduction_sum(double * a, unsigned long size);
*/
        void cuda_scaled_sum_two_float(float * x, const float * y, float b, unsigned long size);
        /*void scaled_sum(float * x, const float * y, const float * z, unsigned long size);
        void scaled_sum(double * x, const double * y, const double * z, unsigned long size);
*/
        void cuda_scale_one_float(const float a, float * x, unsigned long size);

        void cuda_sum_two_float(float * a, const float * b, unsigned long size);
        /*void sum(const float a, float * x, unsigned long size);
        void sum(const double a, double * x, unsigned long size);*/
}
#endif
