
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
        void cuda_difference_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);

        float cuda_dot_product_two_float(const void * a, const void *b, unsigned long size, unsigned long blocksize,
                unsigned long gridsize);

        void cuda_element_inverse_one_float(float * x, unsigned long size, unsigned long blocksize);

        void cuda_element_product_two_float(float * a, const float * b, unsigned long size, unsigned long blocksize);

        void cuda_product_bmdv_q1_float(void * ll, void * ld, void * lu,
                void * dl, void * dd, void *du,
                void * ul, void * ud, void *uu, void * x, void * y,
                unsigned long size, unsigned long blocksize, unsigned long m);

        void cuda_scaled_sum_two_float(void * x, const void * y, float b, unsigned long size, unsigned long blocksize);
        void cuda_scaled_sum_three_float(void * x, void * y, void * z, unsigned long size, unsigned long blocksize);

        void cuda_scale_one_float(void * x, const float a, unsigned long size, unsigned long blocksize);

        void cuda_sum_two_float(void * a, const void * b, unsigned long size, unsigned long blocksize);
}
#endif
