/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBMATH_GUARD_FFT_HH
#define LIBMATH_GUARD_FFT_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_vector.hh>
#include <iostream>
#include "mkl_dfti.h"

namespace honei
{
    struct FFT
    {
        static DenseVector<double> forward(DenseVector<double> & in)
        {
            DenseVector<double> out(in.size() + 2);
            DFTI_DESCRIPTOR_HANDLE my_desc2_handle;
            MKL_LONG status;
            status = DftiCreateDescriptor( &my_desc2_handle, DFTI_DOUBLE,
                    DFTI_REAL, 1, in.size());
            status = DftiSetValue( my_desc2_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
            status = DftiSetValue( my_desc2_handle, DFTI_FORWARD_DOMAIN, DFTI_REAL);
            status = DftiCommitDescriptor( my_desc2_handle);
            status = DftiComputeForward( my_desc2_handle, in.elements(), out.elements());
            status = DftiFreeDescriptor(&my_desc2_handle);

            return out;
        }

        static DenseVector<double> backward(DenseVector<double> & in)
        {
            DenseVector<double> out(in.size() - 2);
            DFTI_DESCRIPTOR_HANDLE my_desc2_handle;
            MKL_LONG status;
            status = DftiCreateDescriptor( &my_desc2_handle, DFTI_DOUBLE,
                    DFTI_REAL, 1, in.size() - 2);
            status = DftiSetValue( my_desc2_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
            status = DftiSetValue( my_desc2_handle, DFTI_FORWARD_DOMAIN, DFTI_REAL);
            status = DftiSetValue( my_desc2_handle, DFTI_BACKWARD_SCALE, 1./(in.size()-2));
            status = DftiCommitDescriptor( my_desc2_handle);
            status = DftiComputeBackward( my_desc2_handle, in.elements(), out.elements());
            status = DftiFreeDescriptor(&my_desc2_handle);

            return out;
        }
    };

}
#endif
