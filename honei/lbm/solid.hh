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

#ifndef LBM_GUARD_SOLID_HH
#define LBM_GUARD_SOLID_HH 1

#include<honei/lbm/tags.hh>
#include<honei/la/dense_matrix.hh>
#include<honei/la/dense_vector.hh>
#include<vector>
#include<algorithm>
#include<cmath>

namespace honei
{
    namespace lbm
    {
        namespace lbm_solid_dims
        {
            class D2;
            class D3;
        }
        namespace lbm_solid_extrapolation_methods
        {
            class SIMPLE;
        }

        template <typename Prec_, typename Dim_>
            class Line
            {
            };

        template <typename Prec_>
            class Line<Prec_, lbm_solid_dims::D2>
            {
                template<typename DT_, typename Dim_>
                    friend class Polygon;

                template<typename Tag_>
                    friend class ScanConversion;

                private:
                Prec_ x_coord_1;
                Prec_ y_coord_1;

                Prec_ x_coord_2;
                Prec_ y_coord_2;

                Prec_ x_min_level;
                Prec_ y_min_level;
                Prec_ x_max_level;
                Prec_ y_max_level;

                public:
                Line<Prec_, lbm_solid_dims::D2>(Prec_ x_1, Prec_ y_1, Prec_ x_2, Prec_ y_2)
                    : x_coord_1(x_1), y_coord_1(y_1), x_coord_2(x_2), y_coord_2(y_2)
                {
                    x_min_level = std::min(x_coord_1, x_coord_2);
                    y_min_level = std::min(y_coord_1, y_coord_2);
                    x_max_level = std::max(x_coord_1, x_coord_2);
                    y_max_level = std::max(y_coord_1, y_coord_2);
                }
            };

        template <typename Prec_, typename Dim_>
            class Polygon
            {
            };

        template<typename Prec_>
            class Polygon<Prec_, lbm_solid_dims::D2>
            {
                template<typename Tag_>
                    friend class ScanConversion;

                private:
                ///Suppose vertices are arranged in line order
                DenseVector<Prec_> * vertex_x_coords;
                DenseVector<Prec_> * vertex_y_coords;
                unsigned long vertex_count, line_count, lines_inserted;
                std::vector< Line<Prec_, lbm_solid_dims::D2> > lines;

                Prec_ line_min_x_level, line_min_y_level, line_max_x_level, line_max_y_level;

                public:
                Polygon<Prec_, lbm_solid_dims::D2>(unsigned long num_lines):
                    vertex_count(num_lines * 2),
                    line_count(num_lines),
                    lines_inserted(0)
                {
                    vertex_x_coords = new DenseVector<Prec_>(vertex_count, Prec_(0));
                    vertex_y_coords = new DenseVector<Prec_>(vertex_count, Prec_(0));
                }

                ~Polygon<Prec_, lbm_solid_dims::D2>()
                {
                    delete vertex_x_coords;
                    delete vertex_y_coords;
                }

                void add_line(Line<Prec_, lbm_solid_dims::D2> & line)
                {
                    if(lines.size() < line_count)
                    {
                        lines.push_back(line);
                        if(lines_inserted == 0)
                        {
                            line_min_x_level = line.x_min_level;
                            line_min_y_level = line.y_min_level;
                            line_max_x_level = line.x_max_level;
                            line_max_y_level = line.y_max_level;
                        }
                        else
                        {
                            line_min_x_level = min(line.x_min_level, line_min_x_level);
                            line_min_y_level = min(line.y_min_level, line_min_y_level);
                            line_max_x_level = max(line.x_max_level, line_max_x_level);
                            line_max_y_level = max(line.y_max_level, line_max_y_level);
                        }
                        ++lines_inserted;
                    }
                    else
                        throw InternalError("Trying to insert " + stringify(lines.size() + 1) + "th line in polygon with " + stringify(line_count) + "lines allowed!");
                }

                void value()
                {
                    unsigned long j(0);
                    for(unsigned long i(0) ; i < line_count * 2; ++i)
                    {
                        (*vertex_x_coords)[i] = lines[j].x_coord_1;
                        (*vertex_y_coords)[i] = lines[j].y_coord_1;
                        ++i;
                        (*vertex_x_coords)[i] = lines[j].x_coord_2;
                        (*vertex_y_coords)[i] = lines[j].y_coord_2;
                        ++j;
                    }
                }

                DenseVector<Prec_> & get_x_coords()
                {
                    return (*vertex_x_coords);
                }

                DenseVector<Prec_> & get_y_coords()
                {
                    return (*vertex_y_coords);
                }

                unsigned long get_line_count()
                {
                    return line_count;
                }
        };


        template<typename Tag_>
            class ScanConversion
            {
            };

        template<>
            class ScanConversion<tags::CPU>
            {
                private:

                    template <typename DT_>
                        static signed long signum(DT_ x)
                        {
                            return (x > 0) ? 1 : (x < 0) ? -1 : 0;
                        }

                    template <typename DT_>
                        static signed long convert_pos(DT_ coord, DT_ delta)
                        {
                            return (signed long)(coord / delta);
                        }

                    ///Bresenham line rasterization:
                    template <typename DT_>
                        static void rasterize_line(Line<DT_, lbm_solid_dims::D2> & line, DenseMatrix<bool> & target, DT_ dx, DT_ dy)
                        {
                            ///Convert to matrix coordinates and determine ranges:
                            signed long d_x(convert_pos(line.x_coord_2, dx) - convert_pos(line.x_coord_1, dx));
                            signed long d_y(convert_pos(line.y_coord_2, dy) - convert_pos(line.y_coord_1, dy));

                            ///Determine signs:
                            signed long inc_x(signum(d_x));
                            signed long inc_y(signum(d_y));
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

                            signed long x(convert_pos(line.x_coord_1, dx));
                            signed long y(convert_pos(line.y_coord_1, dy));
                            signed long err(e_l / 2);

                            ///Set start pixel and begin loop:
                            if(y < target.rows() && x < target.columns() && x >= 0 && y >= 0)
                                target[y][x] = true;

                            for(unsigned long i(0) ; i < e_l ; ++i)
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
                                if(y < target.rows() && x < target.columns() && x >= 0 && y >= 0)
                                    target[y][x] = true;
                            }
                        }

                    ///Local scan fill algo:
                    template <typename DT_>
                        static void local_scan_fill(Polygon<DT_, lbm_solid_dims::D2> & polygon, DenseMatrix<bool> & target, DT_ dx, DT_ dy, bool rect)
                        {
                            signed long i_start_s(convert_pos(polygon.line_min_y_level, dy));
                            signed long j_start_s(convert_pos(polygon.line_min_x_level, dx));
                            signed long i_end_s(convert_pos(polygon.line_max_y_level, dy));
                            signed long j_end_s(convert_pos(polygon.line_max_x_level, dx));

                            unsigned long i_start, i_end, j_start, j_end;
                            if(i_start_s < 0)
                                i_start = 0;
                            else if(i_start_s >= target.rows())
                                i_start = (unsigned long)target.rows();
                            else
                                i_start = (unsigned long)i_start_s;

                            if(i_end_s < 0)
                                i_end = 0;
                            else if(i_end_s >= target.rows())
                                i_end = (unsigned long)target.rows();
                            else
                                i_end = (unsigned long)i_end_s;

                            if(j_start_s < 0)
                                j_start = 0;
                            else if(j_start_s >= target.columns())
                                j_start = (unsigned long)target.columns();
                            else
                                j_start = (unsigned long)j_start_s;

                            if(j_end_s < 0)
                                j_end = 0;
                            else if(j_end_s >= target.columns())
                                j_end = (unsigned long)target.columns();
                            else
                                j_end = (unsigned long)j_end_s;

                            if(rect)
                            {
                                bool a(false), b(false);
                                for(unsigned long i(i_start) ; i < i_end ; ++i)
                                {
                                    for(unsigned long j(j_start); j < j_end ; ++j)
                                    {
                                        bool e(target[i][j]);
                                        //TODO: correct for rectangles, incorrect for triangles
                                        bool rim(j - j_start == 0);
                                        bool a_t((!a & !b & e) | (!a & b & e) | (a & !b & !e) | (a & !b & e) | (a & b & e) | (a & b & !e & !rim));
                                        bool b_t((!a & !b & !e) | (!a & b & !e) | (a & !b & e) | (a & b & !e) | (a & b & e));

                                        target[i][j] = a_t;
                                        a = a_t;
                                        b = b_t;

                                    }
                                }
                            }
                            else
                            {
                                bool a(false), b(false);
                                for(unsigned long i(i_start) ; i < i_end ; ++i)
                                {
                                    for(unsigned long j(j_start); j < j_end ; ++j)
                                    {
                                        bool e(target[i][j]);
                                        bool a_t((!a & !b & e) | (!a & b & e) | (a & !b & !e) | (a & !b & e) | (a & b & e));
                                        bool b_t((!a & !b & !e) | (!a & b & !e) | (a & !b & e) | (a & b & !e) | (a & b & e));

                                        target[i][j] = a_t;
                                        a = a_t;
                                        b = b_t;

                                    }
                                }
                            }
                        }


                public:
                    template<typename Prec_>
                        static void value(Polygon<Prec_, lbm_solid_dims::D2> & solid, DenseMatrix<bool> & target, Prec_ dx, Prec_ dy, bool rect)
                        {
                            ///For all lines: Rasterize line with Bresenhams algo:
                            for(unsigned long i(0) ; i < solid.line_count ; ++i)
                            {
                                rasterize_line(solid.lines[i], target, dx, dy);
                            }

                            ///Fill Polygon:
                            local_scan_fill(solid, target, dx, dy, rect);
                        }
            };

        ///Determiner for fluid cells that become solid
        template<typename Tag_>
            class FluidToSolidCells
            {
            };

        template<>
            class FluidToSolidCells<tags::CPU>
            {
                public:
                    static void value(DenseMatrix<bool> & source_tm1, DenseMatrix<bool> & source_t, DenseMatrix<bool> & target)
                    {
                        for (unsigned long i(0) ; i < target.rows() ; ++i)
                        {
                            for (unsigned long j(0) ; j < target.columns() ; ++j)
                            {
                                target[i][j] = (!source_tm1[i][j] & source_t[i][j]);
                            }
                        }
                    }

            };
        ///Determiner for solid cells that become fluid
        template<typename Tag_>
            class SolidToFluidCells
            {
            };

        template<>
            class SolidToFluidCells<tags::CPU>
            {
                public:
                    static void value(DenseMatrix<bool> & source_tm1, DenseMatrix<bool> & source_t, DenseMatrix<bool> & target)
                    {
                        for (unsigned long i(0) ; i < target.rows() ; ++i)
                        {
                            for (unsigned long j(0) ; j < target.columns() ; ++j)
                            {
                                target[i][j] = (source_tm1[i][j] & !source_t[i][j]);
                            }
                        }
                    }

            };

        ///Initialization of fluids
        template<typename Tag_, typename Method_>
            class STFExtrapolation
            {
            };

        template<>
            class STFExtrapolation<tags::CPU, lbm_solid_extrapolation_methods::SIMPLE>
            {
                public:
                    template<typename Prec_>
                        static void value(DenseMatrix<bool> & obstacles,
                                DenseMatrix<bool> & solids_to_fluids,
                                DenseMatrix<Prec_> & target_h,
                                DenseMatrix<Prec_> & target_u,
                                DenseMatrix<Prec_> & target_v,
                                Prec_ u,     //x-veloc of solid
                                Prec_ v,     //y-veloc of solid
                                Prec_ alpha, //control parameter
                                Prec_ depth) //of solid under water
                        {
                            ///For all cells to init: take average of all fluid neighbours (for h), set vector (for u,v):

                            //First: preprocess elements (only operated on the boolean stf matrix):
                            std::vector<unsigned long> stf_row_index;
                            std::vector<unsigned long> stf_column_index;

                            unsigned long stf_count(0);
                            for(unsigned long i(0) ; i < solids_to_fluids.rows() ; ++i)
                            {
                                for(unsigned long j(0) ; j < solids_to_fluids.columns() ; ++j)
                                {
                                    if(solids_to_fluids[i][j])
                                    {
                                        stf_row_index.push_back(i);
                                        stf_column_index.push_back(j);
                                        ++stf_count;
                                    }
                                }
                            }

                            //Second: determine neighbours and set values:
                            for(unsigned long i(0) ; i < stf_count ; ++i)
                            {
                                std::vector<unsigned long> nb_row_index;
                                std::vector<unsigned long> nb_column_index;

                                unsigned long nb_count(0);
                                //dir 1
                                if(stf_column_index[i] + 1 < obstacles.columns())
                                    if(!obstacles[stf_row_index[i]][stf_column_index[i] + 1] && !solids_to_fluids[stf_row_index[i]][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i]);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //dir 2
                                if(stf_row_index[i] + 1 < obstacles.rows() && stf_column_index[i] + 1 < obstacles.columns())
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i] + 1] && !solids_to_fluids[stf_row_index[i] + 1][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //dir 3
                                if(stf_row_index[i] + 1 < obstacles.rows())
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i]] && !solids_to_fluids[stf_row_index[i] + 1][stf_column_index[i]])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i]);
                                        ++nb_count;
                                    }

                                //dir 4
                                if(stf_row_index[i] + 1 < obstacles.rows() && stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i] - 1] && !solids_to_fluids[stf_row_index[i] + 1][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 5
                                if(stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i]][stf_column_index[i] - 1] && !solids_to_fluids[stf_row_index[i]][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i]);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 6
                                if(stf_row_index[i] - 1 >= 0 && stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i] - 1] && !solids_to_fluids[stf_row_index[i] - 1][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 7
                                if(stf_row_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i]] && !solids_to_fluids[stf_row_index[i] - 1][stf_column_index[i]])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i]);
                                        ++nb_count;
                                    }

                                //dir 8
                                if(stf_row_index[i] - 1 >= 0 && stf_column_index[i] + 1 < obstacles.columns() && !solids_to_fluids[stf_row_index[i] - 1][stf_column_index[i] + 1])
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //determine {h, u, v}_average
                                Prec_ h_average(0);
                                Prec_ u_average(0); //TODO:not needed?
                                Prec_ v_average(0); //TODO:not needed?
                                if(nb_count > 0)
                                {
                                    for(unsigned long j(0) ; j < nb_count ; ++j)
                                    {

                                        Prec_ value_h(alpha * depth * target_h[nb_row_index[j]][nb_column_index[j]]);
                                        Prec_ value_u(alpha * depth * target_u[nb_row_index[j]][nb_column_index[j]]);
                                        Prec_ value_v(alpha * depth * target_v[nb_row_index[j]][nb_column_index[j]]);

                                        h_average += value_h;
                                        u_average += value_u;
                                        v_average += value_v;
                                    }
                                    h_average /= nb_count;
                                    u_average /= nb_count;
                                    v_average /= nb_count;
                                }

                                //determine h, u, v:
                                if(nb_count > 0)
                                {
                                    for(unsigned long j(0) ; j < nb_count ; ++j)
                                    {
                                        target_h[stf_row_index[i]][stf_column_index[i]] = h_average;
                                        target_u[stf_row_index[i]][stf_column_index[i]] = u_average;
                                        target_v[stf_row_index[i]][stf_column_index[i]] = v_average;
                                    }
                                }
                            }
                        }
            };

        ///Initialization of solids
        template<typename Tag_, typename Method_>
            class FTSExtrapolation
            {
            };

        template<>
            class FTSExtrapolation<tags::CPU, lbm_solid_extrapolation_methods::SIMPLE>
            {
                public:
                    template<typename Prec_>
                        static void value(DenseMatrix<bool> & obstacles,
                                DenseMatrix<bool> & fluids_to_solids,
                                DenseMatrix<Prec_> & target_h,
                                DenseMatrix<Prec_> & target_u,
                                DenseMatrix<Prec_> & target_v,
                                Prec_ u,
                                Prec_ v,
                                Prec_ alpha,
                                Prec_ depth)
                        {
                            ///For h: determine displaced amount of water by diving-depth and distribute it amongst FTS cells
                            //First: preprocess elements (only operated on the boolean stf matrix):
                            std::vector<unsigned long> stf_row_index;
                            std::vector<unsigned long> stf_column_index;

                            unsigned long stf_count(0);
                            for(unsigned long i(0) ; i < fluids_to_solids.rows() ; ++i)
                            {
                                for(unsigned long j(0) ; j < fluids_to_solids.columns() ; ++j)
                                {
                                    if(fluids_to_solids[i][j])
                                    {
                                        stf_row_index.push_back(i);
                                        stf_column_index.push_back(j);
                                        ++stf_count;
                                    }
                                }
                            }

                            //Second: determine neighbours and set values:
                            for(unsigned long i(0) ; i < stf_count ; ++i)
                            {
                                std::vector<unsigned long> nb_row_index;
                                std::vector<unsigned long> nb_column_index;

                                unsigned long nb_count(0);
                                //dir 1
                                if(stf_column_index[i] + 1 < obstacles.columns())
                                    if(!obstacles[stf_row_index[i]][stf_column_index[i] + 1] && !fluids_to_solids[stf_row_index[i]][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i]);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //dir 2
                                if(stf_row_index[i] + 1 < obstacles.rows() && stf_column_index[i] + 1 < obstacles.columns())
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i] + 1] && !fluids_to_solids[stf_row_index[i] + 1][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //dir 3
                                if(stf_row_index[i] + 1 < obstacles.rows())
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i]] && !fluids_to_solids[stf_row_index[i] + 1][stf_column_index[i]])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i]);
                                        ++nb_count;
                                    }

                                //dir 4
                                if(stf_row_index[i] + 1 < obstacles.rows() && stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] + 1][stf_column_index[i] - 1] && !fluids_to_solids[stf_row_index[i] + 1][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] + 1);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 5
                                if(stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i]][stf_column_index[i] - 1] && !fluids_to_solids[stf_row_index[i]][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i]);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 6
                                if(stf_row_index[i] - 1 >= 0 && stf_column_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i] - 1] && !fluids_to_solids[stf_row_index[i] - 1][stf_column_index[i] - 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i] - 1);
                                        ++nb_count;
                                    }

                                //dir 7
                                if(stf_row_index[i] - 1 >= 0)
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i]] && !fluids_to_solids[stf_row_index[i] - 1][stf_column_index[i]])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i]);
                                        ++nb_count;
                                    }

                                //dir 8
                                if(stf_row_index[i] - 1 >= 0 && stf_column_index[i] + 1 < obstacles.columns() && !fluids_to_solids[stf_row_index[i] - 1][stf_column_index[i] + 1])
                                    if(!obstacles[stf_row_index[i] - 1][stf_column_index[i] + 1])
                                    {
                                        nb_row_index.push_back(stf_row_index[i] - 1);
                                        nb_column_index.push_back(stf_column_index[i] + 1);
                                        ++nb_count;
                                    }

                                //distribute h, u, v:
                                /*if(nb_count > 0)
                                  {
                                  for(unsigned long j(0) ; j < nb_count ; ++j)
                                  {
                                  target_h[nb_row_index[j]][nb_column_index[j]] += (depth)/nb_count * alpha;
                                  target_u[nb_row_index[j]][nb_column_index[j]] += (depth)/nb_count * alpha;
                                  target_v[nb_row_index[j]][nb_column_index[j]] += (depth)/nb_count * alpha;
                                  }
                                  }*/
                                //TODO: remove these after testing:
                                target_h[stf_row_index[i]][stf_column_index[i]] = Prec_(0);
                                target_u[stf_row_index[i]][stf_column_index[i]] = u;
                                target_v[stf_row_index[i]][stf_column_index[i]] = v;
                            }
                        }
            };

        template<typename Tag_>
            class Positions
            {
            };

        template<>
            class Positions<tags::CPU>
            {
                public:
                    static void value(DenseMatrix<bool> & flags,
                                      std::vector<unsigned long> & row_i,
                                      std::vector<unsigned long> & column_i)
                    {
                        for(unsigned long i(0) ; i < flags.rows() ; ++i)
                        {
                            for(unsigned long j(0) ; j < flags.columns() ; ++j)
                            {
                                if(flags[i][j])
                                {
                                    row_i.push_back(i);
                                    column_i.push_back(j);
                                }
                            }
                        }
                    }
            };

        template<typename Tag_, typename Dir_>
            class Extrapolation
            {
            };

        template<>
            class Extrapolation<tags::CPU, lbm_lattice_types::D2Q9::DIR_1>
            {
                private:
                    template<typename DT_>
                        static DT_ _extrapolation(DenseMatrix<DT_> & target,
                                                  DenseMatrix<bool> & obstacles,
                                                  unsigned long i,
                                                  unsigned long j,
                                                  DT_ dx,
                                                  DT_ h_b)
                        {
                            bool prev((j >= 1) ? true : false);
                            bool pre_prev((j >= 2) ? true : false);

                            bool prev_o(prev ? obstacles[i][j - 1] : true);
                            bool pre_prev_o(pre_prev ? obstacles[i][j - 2] : true);

                            DT_ v_1((prev && !prev_o) ? target[i][j - 1] : h_b);
                            DT_ v_2((pre_prev && !prev_o) ? target[i][j - 2] : (prev && !prev_o) ? target[i][j - 1] : h_b);

                            DT_ m((v_2 - v_1) / dx);
                            DT_ b(v_1 - (m * (j - 1) * dx));

                            return m * j * dx + b;
                        }

                    template<typename DT_>
                        static DT_ _interpolation(DenseMatrix<DT_> & target,
                                                  DenseMatrix<bool> & obstacles,
                                                  unsigned long i,
                                                  unsigned long j,
                                                  DT_ dx,
                                                  DT_ u_x)
                        {
                            bool prev((j >= 1) ? true : false);
                            bool prev_o(prev ? obstacles[i][j - 1] : true);

                            DT_ v_1((prev && !prev_o) ? target[i][j - 1] : DT_(0));

                            DT_ m((u_x - v_1) / DT_(2) * dx);
                            DT_ b(v_1 - (m * (j - 1) * dx));

                            return m * j * dx + b;

                        }
                ///Positive x direction
                public:
                    template<typename DT_>
                        static void value(DenseMatrix<bool> & flags,
                                          DenseMatrix<DT_> & h,
                                          DenseMatrix<DT_> & u,
                                          DenseMatrix<DT_> & v,
                                          DenseMatrix<bool> & obstacles,
                                          DT_ dx,
                                          DT_ u_x,
                                          DT_ h_b)
                        {
                            std::vector<unsigned long> row_i, column_i;

                            Positions<tags::CPU>::value(flags, row_i, column_i);

                            std::vector<unsigned long>::iterator i(row_i.begin());
                            std::vector<unsigned long>::iterator j(column_i.begin());

                            std::vector<unsigned long>::iterator end(row_i.end());

                            for( ; i < end ; ++i, ++j)
                            {
                                h[*i][*j] = _extrapolation(h, obstacles, *i, *j, dx, h_b);
                                u[*i][*j] = _interpolation(u, obstacles, *i, *j, dx, u_x);
                                v[*i][*j] = _interpolation(v, obstacles, *i, *j, dx, DT_(0));
                            }
                        }
            };


    }
}
#endif
