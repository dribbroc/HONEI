/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/spe/liblbm/operations.hh>
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void eq_dist_grid_dir_0_float(vector float * f_eq, const vector float * h, const vector float * u,
                    const vector float * v,
                    const unsigned size, const float g, const float e, const float d_x, const float d_y)
            {
                const vector float g_vector = spu_splats(g);
                const vector float e_vector = spu_splats(e);
                const vector float one = spu_splats(1.f);
                const vector float two = spu_splats(2.f);
                const vector float three = spu_splats(3.f);
                const vector float five = spu_splats(5.f);
                const vector float six = spu_splats(6.f);
                vector float t1, t2, m1, m2, id, td, vd;
                for (unsigned i(0) ; i < size ; ++i)
                {
                    m1 = spu_mul(five, g_vector);
                    m2 = spu_mul(h[i], h[i]);
                    m1 = spu_mul(m1, m2);
                    m2 = spu_mul(e_vector, six);
                    //t1 = divf4(m1, m2);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t1 = spu_mul(m1, vd);

                    m1 = spu_mul(two, h[i]);
                    m2 = spu_mul(e_vector, three);
                    //t2 = divf4(m1, m2);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t2 = spu_mul(m1, vd);

                    m1 = spu_mul(u[i], u[i]);
                    m2 = spu_mul(v[i], v[i]);
                    m1 = spu_add(m1, m2);
                    t2 = spu_mul(t2, m1);

                    m1 = spu_sub(h[i], t1);
                    m1 = spu_sub(m1, t2);
                    f_eq[i] = m1;
                }
            }

            void eq_dist_grid_dir_odd_float(vector float * f_eq, const vector float * h, const vector float * u,
                    const vector float * v,
                    const unsigned size, const float g, const float e, const float d_x, const float d_y)
            {
                const vector float g_vector = spu_splats(g);
                const vector float e_vector = spu_splats(e);
                const vector float d_x_vector = spu_splats(d_x);
                const vector float d_y_vector = spu_splats(d_y);
                const vector float one = spu_splats(1.f);
                const vector float two = spu_splats(2.f);
                const vector float three = spu_splats(3.f);
                const vector float six = spu_splats(6.f);
                vector float t1, t2, t3, t4, dxu, dyv, m1, m2, id, td, vd, scal;
                for (unsigned i(0) ; i < size ; ++i)
                {
                    /*m1 = _mm_set1_ps(distribution_x[dir]);
                    m2 = _mm_load_ps(u + index);
                    dxu = _mm_mul_ps(m1, m2);
                    m1 = _mm_set1_ps(distribution_y[dir]);
                    m2 = _mm_load_ps(v + index);
                    dyv = _mm_mul_ps(m1, m2);*/
                    dxu = spu_mul(d_x_vector, u[i]);
                    dyv = spu_mul(d_y_vector, v[i]);

                    /*scal1 = _mm_set1_ps(float(6));
                    m1 = _mm_load_ps(h + index);
                    m1 = _mm_mul_ps(m1, m1);
                    m1 = _mm_mul_ps(m1, gv);
                    scal1 = _mm_mul_ps(e2v, scal1);
                    t1 = _mm_div_ps(m1, scal1);*/
                    m1 = spu_mul(h[i], h[i]);
                    m1 = spu_mul(m1, g_vector);
                    scal = spu_mul(e_vector, six);
                    vd = scal;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t1 = spu_mul(m1, vd);

                    /*m1 = _mm_load_ps(h + index);
                    t4 = _mm_div_ps(m1, scal1);
                    m1 = _mm_load_ps(u + index);
                    m1 = _mm_mul_ps(m1, m1);
                    m2 = _mm_load_ps(v + index);
                    m2 = _mm_mul_ps(m2, m2);
                    m1 = _mm_add_ps(m1, m2);
                    t4 = _mm_mul_ps(t4, m1);*/
                    vd = scal;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t4 = spu_mul(h[i], vd);
                    m1 = spu_mul(u[i], u[i]);
                    m2 = spu_mul(v[i], v[i]);
                    m1 = spu_add(m1, m2);
                    t4 = spu_mul(t4, m1);

                    /*scal1 = _mm_set1_ps(float(3));
                    m1 = _mm_load_ps(h + index);
                    m2 = _mm_mul_ps(e2v, scal1);
                    m1 = _mm_div_ps(m1, m2);
                    m2 = _mm_add_ps(dxu, dyv);
                    t2 = _mm_mul_ps(m1, m2);*/
                    m2 = spu_mul(e_vector, three);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    m1 = spu_mul(h[i], vd);
                    m2 = spu_add(dxu, dyv);
                    t2 = spu_mul(m1, m2);

                    /*scal1 = _mm_set1_ps(float(2));
                    m1 = _mm_load_ps(h + index);
                    m2 = _mm_mul_ps(e2v, e2v);
                    m2 = _mm_mul_ps(m2, scal1);
                    t3 = _mm_div_ps(m1, m2);*/
                    m2 = spu_mul(e_vector, e_vector);
                    m2 = spu_mul(m2, two);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t3 = spu_mul(h[i], vd);

                    /*m1 = _mm_mul_ps(dxu, dxu);
                    m2 = _mm_mul_ps(scal1, dxu);
                    m2 = _mm_mul_ps(m2, dyv);
                    m1 = _mm_add_ps(m1, m2);
                    m2 = _mm_mul_ps(dyv, dyv);
                    m1 = _mm_add_ps(m1, m2);*/
                    m1 = spu_mul(dxu, dxu);
                    m2 = spu_mul(two, dxu);
                    m2 = spu_mul(m2, dyv);
                    m1 = spu_add(m1, m2);
                    m2 = spu_mul(dyv, dyv);
                    m1 = spu_add(m1, m2);

                    //t3 = _mm_mul_ps(m1, t3);
                    t3 = spu_mul(m1, t3);

                    /*m1 = _mm_add_ps(t1, t2);
                    m1 = _mm_add_ps(m1, t3);
                    m1 = _mm_sub_ps(m1, t4);
                    _mm_store_ps(f_eq + index, m1);*/
                    m1 = spu_add(t1, t2);
                    m1 = spu_add(m1, t3);
                    f_eq[i] = spu_sub(m1, t4);
                }
            }

            void eq_dist_grid_dir_even_float(vector float * f_eq, const vector float * h, const vector float * u,
                    const vector float * v,
                    const unsigned size, const float g, const float e, const float d_x, const float d_y)
            {
                const vector float g_vector = spu_splats(g);
                const vector float e_vector = spu_splats(e);
                const vector float d_x_vector = spu_splats(d_x);
                const vector float d_y_vector = spu_splats(d_y);
                const vector float one = spu_splats(1.f);
                const vector float two = spu_splats(2.f);
                const vector float eight = spu_splats(8.f);
                const vector float twelve = spu_splats(12.f);
                const vector float twentyfour = spu_splats(24.f);
                vector float t1, t2, t3, t4, dxu, dyv, m1, m2, id, td, vd, scal;
                for (unsigned i(0) ; i < size ; ++i)
                {
                    /*m1 = _mm_set1_ps(distribution_x[dir]);
                      m2 = _mm_load_ps(u + index);
                      dxu = _mm_mul_ps(m1, m2);
                      m1 = _mm_set1_ps(distribution_y[dir]);
                      m2 = _mm_load_ps(v + index);
                      dyv = _mm_mul_ps(m1, m2);*/
                    dxu = spu_mul(d_x_vector, u[i]);
                    dyv = spu_mul(d_y_vector, v[i]);

                    /*scal1 = _mm_set1_ps(float(24));
                      m1 = _mm_load_ps(h + index);
                      m1 = _mm_mul_ps(m1, m1);
                      m1 = _mm_mul_ps(m1, gv);
                      scal1 = _mm_mul_ps(e2v, scal1);
                      t1 = _mm_div_ps(m1, scal1);*/
                    m1 = spu_mul(h[i], h[i]);
                    m1 = spu_mul(m1, g_vector);
                    scal = spu_mul(e_vector, twentyfour);
                    vd = scal;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t1 = spu_mul(m1, vd);

                    /*m1 = _mm_load_ps(h + index);
                      t4 = _mm_div_ps(m1, scal1);
                      m1 = _mm_load_ps(u + index);
                      m1 = _mm_mul_ps(m1, m1);
                      m2 = _mm_load_ps(v + index);
                      m2 = _mm_mul_ps(m2, m2);
                      m1 = _mm_add_ps(m1, m2);
                      t4 = _mm_mul_ps(t4, m1);*/
                    vd = scal;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t4 = spu_mul(h[i], vd);
                    m1 = spu_mul(u[i], u[i]);
                    m2 = spu_mul(v[i], v[i]);
                    m1 = spu_add(m1, m2);
                    t4 = spu_mul(t4, m1);

                    /*scal1 = _mm_set1_ps(float(12));
                      m1 = _mm_load_ps(h + index);
                      m2 = _mm_mul_ps(e2v, scal1);
                      m1 = _mm_div_ps(m1, m2);
                      m2 = _mm_add_ps(dxu, dyv);
                      t2 = _mm_mul_ps(m1, m2);*/
                    m2 = spu_mul(e_vector, twelve);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    m1 = spu_mul(h[i], vd);
                    m2 = spu_add(dxu, dyv);
                    t2 = spu_mul(m1, m2);

                    /*scal1 = _mm_set1_ps(float(8));
                      m1 = _mm_load_ps(h + index);
                      m2 = _mm_mul_ps(e2v, e2v);
                      m2 = _mm_mul_ps(m2, scal1);
                      t3 = _mm_div_ps(m1, m2);*/
                    m2 = spu_mul(e_vector, e_vector);
                    m2 = spu_mul(m2, eight);
                    vd = m2;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    t3 = spu_mul(h[i], vd);

                    /*m1 = _mm_mul_ps(dxu, dxu);
                      scal1 = _mm_set1_ps(float(2));
                      m2 = _mm_mul_ps(scal1, dxu);
                      m2 = _mm_mul_ps(m2, dyv);
                      m1 = _mm_add_ps(m1, m2);
                      m2 = _mm_mul_ps(dyv, dyv);
                      m1 = _mm_add_ps(m1, m2);
                      t3 = _mm_mul_ps(m1, t3);*/
                    m1 = spu_mul(dxu, dxu);
                    m2 = spu_mul(dxu, two);
                    m2 = spu_mul(m2, dyv);
                    m1 = spu_add(m1, m2);
                    m2 = spu_mul(dyv, dyv);
                    m1 = spu_add(m1, m2);
                    t3 = spu_mul(m1, t3);

                    /*m1 = _mm_add_ps(t1, t2);
                      m1 = _mm_add_ps(m1, t3);
                      m1 = _mm_sub_ps(m1, t4);
                      _mm_store_ps(f_eq + index, m1);*/
                    m1 = spu_add(t1, t2);
                    m1 = spu_add(m1, t3);
                    f_eq[i] = spu_sub(m1, t4);
                }
            }
        }

        namespace operations
        {
            Operation<4, float, rtm_dma> eq_dist_grid_dir_0_float = {
                &implementation::eq_dist_grid_dir_0_float
            };

            Operation<4, float, rtm_dma> eq_dist_grid_dir_odd_float = {
                &implementation::eq_dist_grid_dir_odd_float
            };

            Operation<4, float, rtm_dma> eq_dist_grid_dir_even_float = {
                &implementation::eq_dist_grid_dir_even_float
            };
        }
    }
}
