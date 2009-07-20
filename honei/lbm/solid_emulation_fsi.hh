/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LBM_GUARD_SOLID_EMULATION_FSI_HH
#define LBM_GUARD_SOLID_EMULATION_FSI_HH 1

#include<honei/lbm/solid.hh>
#include<honei/la/sum.hh>
#include<honei/la/difference.hh>
#include<honei/la/dot_product.hh>
#include<honei/la/product.hh>
#include<honei/la/scale.hh>

namespace honei
{
    namespace lbm
    {
        template<typename Tag_>
            struct MinAngle2D
            {
            };

        template<>
            struct MinAngle2D<tags::CPU>
            {
                public:
                    template<typename DT_>
                        static unsigned long value(DT_ tx, DT_ ty, PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data)

                        {
                            DenseVector<DT_> t(2ul, DT_(0));
                            t[0] = tx;
                            t[0] = ty;
                            if(tx > 0 && ty > 0)
                            {
                                if(tx == ty)
                                    return 6;

                                if(ty < tx)
                                {
                                    ///iterate from 5 to 6
                                    if(ty <= 1./2. * tx)
                                        return 5;
                                    else
                                        return 6;
                                }
                                else
                                {
                                    ///iterate from 6 to 7
                                    if(tx <= 1./2. * ty)
                                        return 7;
                                    else
                                        return 6;
                                }
                            }
                            if(tx < 0 && ty < 0)
                            {
                                if(tx == ty)
                                    return 2;

                                if(ty < tx)
                                {
                                    if(ty <= 1./2. * tx)
                                        return 3;
                                    else
                                        return 2;
                                }
                                else
                                {
                                    if(tx <= 1./2. * ty)
                                        return 1;
                                    else
                                        return 2;
                                }
                            }
                            if(tx > 0 && ty < 0)
                            {
                                if(tx == ty)
                                    return 4;

                                if(-ty < tx)
                                {
                                    ///iterate from 4 to 5
                                    if(-ty <= 1./2. * tx)
                                        return 5;
                                    else
                                        return 4;
                                }
                                else
                                {
                                    ///iterate from 3 to 4
                                    if(tx <= 1./2. * -ty)
                                        return 3;
                                    else
                                        return 4;
                                }
                            }
                            if(tx < 0 && ty > 0)
                            {
                                if(tx == -ty)
                                    return 8;

                                if(-tx < ty)
                                {
                                    ///iterate from 8 to 7
                                    if(-tx <= 1./2. * ty)
                                        return 7;
                                    else
                                        return 8;
                                }
                                else
                                {
                                    ///iterate from 8 to 1
                                    if(-tx >= 1./2. * ty)
                                        return 1;
                                    else
                                        return 8;
                                }
                            }

                            if(tx == 0)
                            {
                                if(ty > 0)
                                    return 7;
                                else
                                    return 3;
                            }
                            if(ty == 0)
                            {
                                if(tx > 0)
                                    return 5;
                                else
                                    return 1;
                            }
                        }
            };

        template<typename Tag_>
            struct SolidEmulation2D
            {
            };

        template<>
            struct SolidEmulation2D<tags::CPU>
            {
                private:
                    template<typename DT_>
                        static DenseVector<DT_> _translation_vector(Polygon<DT_, lbm_solid_dims::D2> & polygon,
                                PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                Grid<D2Q9, DT_> & grid)
                        {
                            ///Accumulate f_mea
                            DenseVector<DT_> F(solids.f_mea_1->size(), DT_(0));
                            Sum<tags::CPU>::value(F, *solids.f_mea_1);
                            Sum<tags::CPU>::value(F, *solids.f_mea_2);
                            Sum<tags::CPU>::value(F, *solids.f_mea_3);
                            Sum<tags::CPU>::value(F, *solids.f_mea_4);
                            Sum<tags::CPU>::value(F, *solids.f_mea_5);
                            Sum<tags::CPU>::value(F, *solids.f_mea_6);
                            Sum<tags::CPU>::value(F, *solids.f_mea_7);
                            Sum<tags::CPU>::value(F, *solids.f_mea_8);

                            solids.f_mea_1->lock(lm_read_and_write);
                            solids.f_mea_2->lock(lm_read_and_write);
                            solids.f_mea_3->lock(lm_read_and_write);
                            solids.f_mea_4->lock(lm_read_and_write);
                            solids.f_mea_5->lock(lm_read_and_write);
                            solids.f_mea_6->lock(lm_read_and_write);
                            solids.f_mea_7->lock(lm_read_and_write);
                            solids.f_mea_8->lock(lm_read_and_write);

                            solids.line_flags->lock(lm_read_only);

                            for(unsigned long i(0) ;  i < F.size() ; ++i)
                            {
                                if(!(*solids.line_flags)[i])
                                    F[i] = DT_(0);
                            }

                            solids.line_flags->unlock(lm_read_only);

                            solids.f_mea_1->unlock(lm_read_and_write);
                            solids.f_mea_2->unlock(lm_read_and_write);
                            solids.f_mea_3->unlock(lm_read_and_write);
                            solids.f_mea_4->unlock(lm_read_and_write);
                            solids.f_mea_5->unlock(lm_read_and_write);
                            solids.f_mea_6->unlock(lm_read_and_write);
                            solids.f_mea_7->unlock(lm_read_and_write);
                            solids.f_mea_8->unlock(lm_read_and_write);
                            DenseMatrix<bool> line_res(grid.h->rows(), grid.h->columns());
                            GridPackerFSI<D2Q9, NOSLIP, DT_>::deflate(grid, solids, solids.line_flags, &line_res);

                            DenseVector<DT_> result(4ul, DT_(0));
                            DenseVector<DT_> point(4ul, DT_(0));
                            DenseVector<DT_> mass(4ul, DT_(0));
                            mass[0] = polygon._mass_x * grid.d_x;
                            mass[1] = polygon._mass_y * grid.d_y;
                            mass[2] = DT_(0);
                            mass[3] = DT_(1);
                            for(unsigned long i(0) ; i < line_res.rows() ; ++i)
                            {
                                for(unsigned long j(0) ; j < line_res.columns() ; ++j)
                                {
                                    if(line_res[i][j])
                                    {
                                        point[0] = j * grid.d_x;
                                        point[1] = i * grid.d_y;
                                        point[2] = DT_(0);
                                        point[3] = DT_(1);

                                        Difference<tags::CPU>::value(point, mass);
                                        Scale<tags::CPU>::value(point, DT_(-1 * F[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, i, j)]));
                                        Sum<tags::CPU>::value(result, point);
                                    }
                                }
                            }
                            return result;

                        }

                    template<typename DT_>
                        static DenseMatrix<DT_> _translation_matrix(DenseVector<DT_> & tv)
                        {
                            DenseMatrix<DT_> result(4ul, 4ul, DT_(0));
                            result[0][0] = DT_(1);
                            result[1][1] = DT_(1);
                            result[2][2] = DT_(1);
                            result[3][3] = DT_(1);

                            result[0][3] = tv[0];
                            result[1][3] = tv[1];
                            result[2][3] = tv[2];

                            return result;
                        }

                public:
                    template<typename DT_>
                        static void value(Polygon<DT_, lbm_solid_dims::D2> & polygon,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          Grid<D2Q9, DT_> & grid)

                        {
                            DenseVector<DT_> tv(_translation_vector(polygon, solids, grid));
                            unsigned long m(MinAngle2D<tags::CPU>::value(tv[0], tv[1], solids, data));
                            tv[0] = (*data.distribution_x)[m];
                            tv[1] = (*data.distribution_y)[m];

                            solids.current_u = grid.d_x * tv[0];
                            solids.current_v = grid.d_y * tv[1];
                            polygon._mass_x = polygon._mass_x + tv[0];
                            polygon._mass_y = polygon._mass_y + tv[1];
                            for(unsigned long i(0); i < polygon.line_count; ++i)
                            {
                                polygon.lines[i].x_coord_1 = polygon.lines[i].x_coord_1 + tv[0];
                                polygon.lines[i].y_coord_1 = polygon.lines[i].y_coord_1 + tv[1];
                                polygon.lines[i].x_coord_2 = polygon.lines[i].x_coord_2 + tv[0];
                                polygon.lines[i].y_coord_2 = polygon.lines[i].y_coord_2 + tv[1];

                                polygon.lines[i].x_min_level = polygon.lines[i].x_min_level + tv[0];
                                polygon.lines[i].y_min_level = polygon.lines[i].y_min_level + tv[1];
                            }
                        }
            };
    }
}
#endif
