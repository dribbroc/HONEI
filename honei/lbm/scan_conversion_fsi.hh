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
#ifndef LBM_GUARD_SCAN_CONVERSION_FSI_HH
#define LBM_GUARD_SCAN_CONVERSION_FSI_HH 1

#include <honei/lbm/tags.hh>
#include <honei/la/algorithm.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/solid.hh>

namespace honei
{
    namespace lbm
    {

        template<typename Tag_>
            class ScanConversionFSI
            {
            };

        template<>
            class ScanConversionFSI<tags::CPU>
            {
                private:

                    template <typename DT_>
                        static inline signed long _signum(DT_ x)
                        {
                            return (x > 0) ? 1 : (x < 0) ? -1 : 0;
                        }

                    template <typename DT_>
                        static inline signed long _convert_pos(DT_ coord, DT_ delta)
                        {
                            return (signed long)(coord / delta);
                        }

                    template <typename DT_>
                        static inline void _flag_line_neighbours(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                                                 PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                                                 PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                                                 unsigned long packed_index)
                        {
                            unsigned long nb_1(((*info.cuda_dir_1)[packed_index] != 4294967295) ? (*info.cuda_dir_1)[packed_index] : packed_index);
                            unsigned long nb_2(((*info.cuda_dir_2)[packed_index] != 4294967295) ? (*info.cuda_dir_2)[packed_index] : packed_index);
                            unsigned long nb_3(((*info.cuda_dir_3)[packed_index] != 4294967295) ? (*info.cuda_dir_3)[packed_index] : packed_index);
                            unsigned long nb_4(((*info.cuda_dir_4)[packed_index] != 4294967295) ? (*info.cuda_dir_4)[packed_index] : packed_index);
                            unsigned long nb_5(((*info.cuda_dir_5)[packed_index] != 4294967295) ? (*info.cuda_dir_5)[packed_index] : packed_index);
                            unsigned long nb_6(((*info.cuda_dir_6)[packed_index] != 4294967295) ? (*info.cuda_dir_6)[packed_index] : packed_index);
                            unsigned long nb_7(((*info.cuda_dir_7)[packed_index] != 4294967295) ? (*info.cuda_dir_7)[packed_index] : packed_index);
                            unsigned long nb_8(((*info.cuda_dir_8)[packed_index] != 4294967295) ? (*info.cuda_dir_8)[packed_index] : packed_index);

                            (*solids.boundary_flags)[nb_1] = true;
                            (*solids.boundary_flags)[nb_2] = true;
                            (*solids.boundary_flags)[nb_3] = true;
                            (*solids.boundary_flags)[nb_4] = true;
                            (*solids.boundary_flags)[nb_5] = true;
                            (*solids.boundary_flags)[nb_6] = true;
                            (*solids.boundary_flags)[nb_7] = true;
                            (*solids.boundary_flags)[nb_8] = true;
                        }

                    template <typename DT_>
                    static void _clamp(Grid<D2Q9, DT_> & grid,
                                       PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                       PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                       PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                       signed long i,
                                       signed long j)
                    {
                        bool north((i < 0) ? true : false);
                        bool west((j < 0) ? true : false);
                        bool south((i >= (signed long)grid.h->rows()) ? true : false);
                        bool east((j >= (signed long)grid.h->columns()) ? true : false);

                        unsigned long target_x(west ? 0ul : east ? grid.h->columns() - 1 : (unsigned long)j);
                        unsigned long target_y(north ? 0ul : south ? grid.h->rows() - 1 : (unsigned long)i);

                        unsigned long packed_index(GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, target_y, target_x));
                        _flag_line_neighbours(info, data, solids, packed_index);
                    }


                    template <typename DT_>
                        static void _rasterize_line(Grid<D2Q9, DT_> & grid,
                                                    PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                                    PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                                    PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                                    Line<DT_, lbm_solid_dims::D2> & line)
                        {
                            DT_ dx(grid.d_x);
                            DT_ dy(grid.d_y);

                            ///Convert to matrix coordinates and determine ranges:
                            signed long d_x(_convert_pos(line.x_coord_2, dx) - _convert_pos(line.x_coord_1, dx));
                            signed long d_y(_convert_pos(line.y_coord_2, dy) - _convert_pos(line.y_coord_1, dy));
                            ///Determine signs:
                            signed long inc_x(_signum(d_x));
                            signed long inc_y(_signum(d_y));
                            if(d_x < 0) d_x = -d_x;
                            if(d_y < 0) d_y = -d_y;

                            signed long p_d_x, p_d_y, d_d_x, d_d_y, e_s, e_l;
                            if(d_x > d_y)
                            {
                                p_d_x = inc_x;
                                p_d_y = 0;
                                d_d_x = inc_x;
                                d_d_y = inc_y;
                                e_s = d_y;
                                e_l = d_x;
                            }
                            else
                            {
                                p_d_x = 0;
                                p_d_y = inc_y;
                                d_d_x = inc_x;
                                d_d_y = inc_y;
                                e_s = d_x;
                                e_l = d_y;
                            }

                            signed long x(_convert_pos(line.x_coord_1, dx));
                            signed long y(_convert_pos(line.y_coord_1, dy));
                            signed long err(e_l / 2);

                            ///Set start pixel and begin loop:
                            if(y < (signed long)grid.h->rows() && x < (signed long)grid.h->columns() && x >= 0 && y >= 0)
                            {
                                (*solids.line_flags)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, y, x)] = true;
                                _flag_line_neighbours(info, data, solids, GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, y, x));
                            }

                            for(signed long i(0) ; i < e_l ; ++i)
                            {
                                err -= e_s;
                                if(err < 0)
                                {
                                    err += e_l;
                                    x += d_d_x;
                                    y += d_d_y;
                                }
                                else
                                {
                                    x += p_d_x;
                                    y += p_d_y;
                                }
                                if(y < (signed long)grid.h->rows() && x < (signed long)grid.h->columns() && x >= 0 && y >= 0)
                                {
                                    (*solids.line_flags)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, y, x)] = true;
                                    _flag_line_neighbours(info, data, solids, GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, y, x));
                                }
                                else
                                {
                                    _clamp(grid, info, data, solids, y, x);
                                }
                            }
                        }

                    template <typename DT_>
                        static void _local_scan_fill(Grid<D2Q9, DT_> & grid,
                                                     PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                                     PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                                     PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                                     Polygon<DT_, lbm_solid_dims::D2> & polygon,
                                                     bool rect)
                        {
                            DT_ dx(grid.d_x);
                            DT_ dy(grid.d_y);

                            signed long i_start_s(_convert_pos(polygon.line_min_y_level, dy));
                            signed long j_start_s(_convert_pos(polygon.line_min_x_level, dx));
                            signed long i_end_s(_convert_pos(polygon.line_max_y_level, dy));
                            signed long j_end_s(_convert_pos(polygon.line_max_x_level, dx));

                            unsigned long i_start, i_end, j_start, j_end;
                            if(i_start_s < 0)
                                i_start = 0;
                            else if(i_start_s >= (signed long)grid.h->rows())
                                i_start = (unsigned long)grid.h->rows();
                            else
                                i_start = (unsigned long)i_start_s;

                            if(i_end_s < 0)
                                i_end = 0;
                            else if(i_end_s >= (signed long)grid.h->rows())
                                i_end = (unsigned long)grid.h->rows();
                            else
                                i_end = (unsigned long)i_end_s;

                            if(j_start_s < 0)
                                j_start = 0;
                            else if(j_start_s >= (signed long)grid.h->columns())
                                j_start = (unsigned long)grid.h->columns();
                            else
                                j_start = (unsigned long)j_start_s;

                            if(j_end_s < 0)
                                j_end = 0;
                            else if(j_end_s >= (signed long)grid.h->columns())
                                j_end = (unsigned long)grid.h->columns();
                            else
                                j_end = (unsigned long)j_end_s;

                            bool a(false), b(false);
                            for(unsigned long i(i_start) ; i <= i_end ; ++i)
                            {
                                for(unsigned long j(j_start); j <= j_end ; ++j)
                                {
                                    bool e_1((*solids.line_flags)[GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, i, j)]);
                                    bool e_2(j - j_start == 1);

                                    bool a_t( (a & b & !e_1 & !e_2) |
                                             (!a & !b & e_1 & !e_2) |
                                             (a & !b & e_1 & !e_2) |
                                             (!a & b & e_1 & !e_2) |
                                             (a & b & e_1 & !e_2) |
                                             (a & !b & !e_1 & e_2) |
                                             (a & b & !e_1 & e_2) |
                                             (!a & !b & e_1 & e_2) |
                                             (a & !b & e_1 & e_2) |
                                             (!a & b & e_1 & e_2) |
                                             (a & b & e_1 & e_2)
                                           );

                                    bool b_t( (rect) ? (!a & !b & !e_1 & !e_2) |
                                              (a & !b & !e_1 & !e_2) |
                                              (!a & b & !e_1 & !e_2) |
                                              (a & b & !e_1 & !e_2) |
                                              (!a & !b & !e_1 & e_2) |
                                              (a & !b & !e_1 & e_2) |
                                              (!a & b & !e_1 & e_2) |
                                              (a & b & !e_1 & e_2) |
                                              (!a & b & e_1 & !e_2) |
                                              (!a & b & e_1 & e_2) :

                                              (j == j_start || j == j_end || i == i_end || i == i_start) ?
                                              (!a & !b & !e_1 & !e_2) |
                                              (a & !b & !e_1 & !e_2) |
                                              (!a & b & !e_1 & !e_2) |
                                              (a & b & !e_1 & !e_2) |
                                              (!a & !b & !e_1 & e_2) |
                                              (a & !b & !e_1 & e_2) |
                                              (!a & b & !e_1 & e_2) |
                                              (a & b & !e_1 & e_2) :

                                              (!a & !b & !e_1 & !e_2) |
                                              (a & !b & !e_1 & !e_2) |
                                              (!a & b & !e_1 & !e_2) |
                                              (a & b & !e_1 & !e_2) |
                                              (!a & !b & !e_1 & e_2) |
                                              (a & !b & !e_1 & e_2) |
                                              (!a & b & !e_1 & e_2) |
                                              (a & b & !e_1 & e_2) |
                                              (!a & b & e_1 & !e_2) |
                                              (!a & b & e_1 & e_2)
                                            );

                                    unsigned long packed_index(GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, i, j));
                                    (*solids.solid_flags)[packed_index] = a_t;
                                    (*solids.boundary_flags)[packed_index] = a_t ? false : (*solids.boundary_flags)[packed_index];

                                    a = a_t;
                                    b = b_t;

                                }
                            }
                        }
                public:
                    template<typename DT_>
                        static void value(Grid<D2Q9, DT_> & grid,
                                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                                          PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
                                          PackedSolidData<lbm_lattice_types::D2Q9, DT_> & solids,
                                          Polygon<DT_, lbm_solid_dims::D2> & solid,
                                          bool rect)
                        {
                            info.cuda_dir_1->lock(lm_read_only);
                            info.cuda_dir_2->lock(lm_read_only);
                            info.cuda_dir_3->lock(lm_read_only);
                            info.cuda_dir_4->lock(lm_read_only);
                            info.cuda_dir_5->lock(lm_read_only);
                            info.cuda_dir_6->lock(lm_read_only);
                            info.cuda_dir_7->lock(lm_read_only);
                            info.cuda_dir_8->lock(lm_read_only);

                            ///Clear all target vectors:
                            fill<tags::CPU>((*solids.line_flags));
                            fill<tags::CPU>((*solids.boundary_flags));
                            fill<tags::CPU>((*solids.solid_flags));
                            solids.line_flags->lock(lm_read_and_write);
                            solids.boundary_flags->lock(lm_read_and_write);
                            solids.solid_flags->lock(lm_read_and_write);

                            grid.h->lock(lm_read_only);


                            ///For all lines: Rasterize line with Bresenhams algo:
                            for(unsigned long i(0) ; i < solid.line_count ; ++i)
                            {
                                _rasterize_line(grid, info, data, solids, solid.lines[i]);
                            }

                            ///Fill Polygon:
                            _local_scan_fill(grid, info, data, solids, solid, rect);

                            solids.line_flags->unlock(lm_read_and_write);
                            solids.boundary_flags->unlock(lm_read_and_write);
                            solids.solid_flags->unlock(lm_read_and_write);

                            info.cuda_dir_1->unlock(lm_read_only);
                            info.cuda_dir_2->unlock(lm_read_only);
                            info.cuda_dir_3->unlock(lm_read_only);
                            info.cuda_dir_4->unlock(lm_read_only);
                            info.cuda_dir_5->unlock(lm_read_only);
                            info.cuda_dir_6->unlock(lm_read_only);
                            info.cuda_dir_7->unlock(lm_read_only);
                            info.cuda_dir_8->unlock(lm_read_only);


                            grid.h->unlock(lm_read_only);
                        }
            };
    }
}

#endif
