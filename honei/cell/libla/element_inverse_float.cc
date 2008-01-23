/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/cell/cell.hh>
#include <honei/cell/libla/operations.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void element_inverse_float(vector float * elements, const unsigned size, const float scalar)
            {
                const Vector<float> one = { { 1.0f, 1.0f, 1.0f, 1.0f } };
                const Vector<float> zero = { { 0.0f, 0.0f, 0.0f, 0.0f } };
                // Bitmask for the exponent field in single precision IEEE 754 floats.
                const Vector<unsigned> mask = { { 0x7f800000U, 0x7f800000U, 0x7f800000U, 0x7f800000U } };

                for (unsigned k(0) ; k < size ; ++k)
                {
                    // Do FREST, FI and one step of Newton-Raphson as proposed
                    // in SPUISAv1.2, p. 215f.
                    Vector<float> v = { elements[k] };
                    Vector<float> i = { spu_re(v.vf) };
                    Vector<float> t = { spu_nmsub(v.vf, i.vf, one.vf) };
                    v.vf = spu_madd(t.vf, i.vf, i.vf);

                    // Replace NANs by zeros.
                    vector unsigned int finite_test(spu_and(v.vui, mask.vui));
                    finite_test = spu_cmpeq(finite_test, mask.vui);
                    v.vui = spu_sel(v.vui, zero.vui, finite_test);
                    elements[k] = v.vf;
                }
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_dma> element_inverse_float = {
                &implementation::element_inverse_float,
            };
        }
    }
}

