/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LibLBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBLBM_GUARD_SOLVER_LABSWE_HH
#define LIBLBM_GUARD_SOLVER_LABSWE_HH 1

/**
 * \file
 * Implementation of a SWE solver using LBM.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sum.hh>
#include <honei/la/scale.hh>
#include <honei/la/product.hh>
#include <honei/la/element_inverse.hh>
#include <cmath>

using namespace honei::lbm;

namespace honei
{

    template<typename Tag_,
        typename ResPrec_,
        typename SourceType_,
        typename SourceSheme_,
        typename GridType_,
        typename LatticeType_,
        typename BoundaryType_>
            class SolverLABSWE
            {
            };

    template<typename Tag_, typename ResPrec_>
        class SolverLABSWE<Tag_, ResPrec_, lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC>
        {
            private:
                /** Global variables.
                 *
                 **/
                ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t;
                DenseMatrix<ResPrec_>* _height;
                DenseMatrix<ResPrec_>* _bottom;
                DenseMatrix<ResPrec_>* _d_bottom_x;
                DenseMatrix<ResPrec_>* _d_bottom_y;
                DenseMatrix<ResPrec_>* _u;
                DenseMatrix<ResPrec_>* _v;

                DenseMatrix<ResPrec_>* _distribution_0;
                DenseMatrix<ResPrec_>* _distribution_1;
                DenseMatrix<ResPrec_>* _distribution_2;
                DenseMatrix<ResPrec_>* _distribution_3;
                DenseMatrix<ResPrec_>* _distribution_4;
                DenseMatrix<ResPrec_>* _distribution_5;
                DenseMatrix<ResPrec_>* _distribution_6;
                DenseMatrix<ResPrec_>* _distribution_7;
                DenseMatrix<ResPrec_>* _distribution_8;

                DenseMatrix<ResPrec_>* _temp_distribution_0;
                DenseMatrix<ResPrec_>* _temp_distribution_1;
                DenseMatrix<ResPrec_>* _temp_distribution_2;
                DenseMatrix<ResPrec_>* _temp_distribution_3;
                DenseMatrix<ResPrec_>* _temp_distribution_4;
                DenseMatrix<ResPrec_>* _temp_distribution_5;
                DenseMatrix<ResPrec_>* _temp_distribution_6;
                DenseMatrix<ResPrec_>* _temp_distribution_7;
                DenseMatrix<ResPrec_>* _temp_distribution_8;

                DenseMatrix<ResPrec_>* _eq_distribution_0;
                DenseMatrix<ResPrec_>* _eq_distribution_1;
                DenseMatrix<ResPrec_>* _eq_distribution_2;
                DenseMatrix<ResPrec_>* _eq_distribution_3;
                DenseMatrix<ResPrec_>* _eq_distribution_4;
                DenseMatrix<ResPrec_>* _eq_distribution_5;
                DenseMatrix<ResPrec_>* _eq_distribution_6;
                DenseMatrix<ResPrec_>* _eq_distribution_7;
                DenseMatrix<ResPrec_>* _eq_distribution_8;

                unsigned long _grid_width, _grid_height;

                DenseMatrix<ResPrec_>* _source_x;
                DenseMatrix<ResPrec_>* _source_y;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi;
                DenseVector<ResPrec_>* _distribution_vector_x;
                DenseVector<ResPrec_>* _distribution_vector_y;

                /** Capsule for the extration of SWE physical quantities.
                 *
                 **/
                void _extract()
                {
                    DenseMatrix<ResPrec_> accu(_distribution_0->copy());

                    Sum<Tag_>::value(accu, *_distribution_1);
                    Sum<Tag_>::value(accu, *_distribution_2);
                    Sum<Tag_>::value(accu, *_distribution_3);
                    Sum<Tag_>::value(accu, *_distribution_4);
                    Sum<Tag_>::value(accu, *_distribution_5);
                    Sum<Tag_>::value(accu, *_distribution_6);
                    Sum<Tag_>::value(accu, *_distribution_7);
                    Sum<Tag_>::value(accu, *_distribution_8);

                    *_height = accu;

                    DenseMatrix<ResPrec_> d0c(_distribution_0->copy());
                    DenseMatrix<ResPrec_> d1c(_distribution_1->copy());
                    DenseMatrix<ResPrec_> d2c(_distribution_2->copy());
                    DenseMatrix<ResPrec_> d3c(_distribution_3->copy());
                    DenseMatrix<ResPrec_> d4c(_distribution_4->copy());
                    DenseMatrix<ResPrec_> d5c(_distribution_5->copy());
                    DenseMatrix<ResPrec_> d6c(_distribution_6->copy());
                    DenseMatrix<ResPrec_> d7c(_distribution_7->copy());
                    DenseMatrix<ResPrec_> d8c(_distribution_8->copy());

                    Scale<Tag_>::value( *_distribution_0, (*_distribution_vector_x)[0]);
                    Scale<Tag_>::value( *_distribution_1, (*_distribution_vector_x)[1]);
                    Scale<Tag_>::value( *_distribution_2, (*_distribution_vector_x)[2]);
                    Scale<Tag_>::value( *_distribution_3, (*_distribution_vector_x)[3]);
                    Scale<Tag_>::value( *_distribution_4, (*_distribution_vector_x)[4]);
                    Scale<Tag_>::value( *_distribution_5, (*_distribution_vector_x)[5]);
                    Scale<Tag_>::value( *_distribution_6, (*_distribution_vector_x)[6]);
                    Scale<Tag_>::value( *_distribution_7, (*_distribution_vector_x)[7]);
                    Scale<Tag_>::value( *_distribution_8, (*_distribution_vector_x)[8]);

                    DenseMatrix<ResPrec_> accu2(_distribution_0->copy());

                    Sum<Tag_>::value(accu2, *_distribution_1);
                    Sum<Tag_>::value(accu2, *_distribution_2);
                    Sum<Tag_>::value(accu2, *_distribution_3);
                    Sum<Tag_>::value(accu2, *_distribution_4);
                    Sum<Tag_>::value(accu2, *_distribution_5);
                    Sum<Tag_>::value(accu2, *_distribution_6);
                    Sum<Tag_>::value(accu2, *_distribution_7);
                    Sum<Tag_>::value(accu2, *_distribution_8);

                    DenseMatrix<ResPrec_> h_inv(_height->copy());
                    ElementInverse<Tag_>::value(h_inv);
                    *_u = Product<Tag_>::value(h_inv, accu2);

                    Scale<Tag_>::value( d0c, (*_distribution_vector_y)[0]);
                    Scale<Tag_>::value( d1c, (*_distribution_vector_y)[1]);
                    Scale<Tag_>::value( d2c, (*_distribution_vector_y)[2]);
                    Scale<Tag_>::value( d3c, (*_distribution_vector_y)[3]);
                    Scale<Tag_>::value( d4c, (*_distribution_vector_y)[4]);
                    Scale<Tag_>::value( d5c, (*_distribution_vector_y)[5]);
                    Scale<Tag_>::value( d6c, (*_distribution_vector_y)[6]);
                    Scale<Tag_>::value( d7c, (*_distribution_vector_y)[7]);
                    Scale<Tag_>::value( d8c, (*_distribution_vector_y)[8]);

                    DenseMatrix<ResPrec_> accu3(d0c.copy());

                    Sum<Tag_>::value(accu3, d1c);
                    Sum<Tag_>::value(accu3, d2c);
                    Sum<Tag_>::value(accu3, d3c);
                    Sum<Tag_>::value(accu3, d4c);
                    Sum<Tag_>::value(accu3, d5c);
                    Sum<Tag_>::value(accu3, d6c);
                    Sum<Tag_>::value(accu3, d7c);
                    Sum<Tag_>::value(accu3, d8c);

                    ElementInverse<Tag_>::value(h_inv);
                    *_v = Product<Tag_>::value(h_inv, accu3);

                }

            public:
                SolverLABSWE(ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, unsigned long gx, unsigned long gy, DenseMatrix<ResPrec_>* height,
                                                                                                        DenseMatrix<ResPrec_>* bottom,
                                                                                                        DenseMatrix<ResPrec_>* u,
                                                                                                        DenseMatrix<ResPrec_>* v) :
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _grid_width(gx),
                    _grid_height(gy),
                    _height(height),
                    _bottom(bottom),
                    _u(u),
                    _v(v),
                    _pi(3.14159265),
                    _gravity(9.81),
                    _n_alpha(ResPrec_(6))
            {
                CONTEXT("When creating LABSWE solver:");
                _e = _delta_x / _delta_t;
                _distribution_vector_x = new DenseVector<ResPrec_>(9ul);
                _distribution_vector_y = new DenseVector<ResPrec_>(9ul);
                _d_bottom_x = new DenseMatrix<ResPrec_>(gx, gy);
                _d_bottom_y = new DenseMatrix<ResPrec_>(gx, gy);

                _distribution_0 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_1 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_2 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_3 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_4 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_5 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_6 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_7 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _distribution_8 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));

                _temp_distribution_0 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_1 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_2 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_3 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_4 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_5 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_6 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_7 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _temp_distribution_8 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));

                _eq_distribution_0 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_1 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_2 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_3 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_4 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_5 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_6 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_7 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));
                _eq_distribution_8 = new DenseMatrix<ResPrec_> (_grid_height, _grid_width, ResPrec_(0));

                (*_distribution_vector_x)[0] = ResPrec_(0.);
                (*_distribution_vector_x)[1] = ResPrec_(_e * cos(0.));
                (*_distribution_vector_x)[2] = ResPrec_(sqrt(2.) * _e * cos(_pi / 4.));
                (*_distribution_vector_x)[3] = ResPrec_(_e * cos(_pi / 2.));
                (*_distribution_vector_x)[4] = ResPrec_(sqrt(2.) * _e * cos(3. * _pi / 4.));
                (*_distribution_vector_x)[5] = ResPrec_(_e * cos(_pi));
                (*_distribution_vector_x)[6] = ResPrec_(sqrt(2.) * _e * cos(5. * _pi / 4.));
                (*_distribution_vector_x)[7] = ResPrec_(_e * cos(3. * _pi / 2.));
                (*_distribution_vector_x)[8] = ResPrec_(sqrt(2.) * _e * cos(7. * _pi / 4.));
                (*_distribution_vector_y)[0] = ResPrec_(0.);
                (*_distribution_vector_y)[1] = ResPrec_(_e * sin(0.));
                (*_distribution_vector_y)[2] = ResPrec_(sqrt(2.) * _e * sin(_pi / 4.));
                (*_distribution_vector_y)[3] = ResPrec_(_e * sin(_pi / 2.));
                (*_distribution_vector_y)[4] = ResPrec_(sqrt(2.) * _e * sin(3. * _pi / 4.));
                (*_distribution_vector_y)[5] = ResPrec_(_e * sin(_pi));
                (*_distribution_vector_y)[6] = ResPrec_(sqrt(2.) * _e * sin(5. * _pi / 4.));
                (*_distribution_vector_y)[7] = ResPrec_(_e * sin(3. * _pi / 2.));
                (*_distribution_vector_y)[8] = ResPrec_(sqrt(2.) * _e * sin(7. * _pi / 4.));
            }

            ~SolverLABSWE()
            {
                CONTEXT("When destroying LABSWE solver.");
                delete _distribution_vector_x;
                delete _distribution_vector_y;
                delete _d_bottom_x;
                delete _d_bottom_y;

                delete _distribution_0;
                delete _distribution_1;
                delete _distribution_2;
                delete _distribution_3;
                delete _distribution_4;
                delete _distribution_5;
                delete _distribution_6;
                delete _distribution_7;
                delete _distribution_8;
                delete _temp_distribution_0;
                delete _temp_distribution_1;
                delete _temp_distribution_2;
                delete _temp_distribution_3;
                delete _temp_distribution_4;
                delete _temp_distribution_5;
                delete _temp_distribution_6;
                delete _temp_distribution_7;
                delete _temp_distribution_8;
                delete _eq_distribution_0;
                delete _eq_distribution_1;
                delete _eq_distribution_2;
                delete _eq_distribution_3;
                delete _eq_distribution_4;
                delete _eq_distribution_5;
                delete _eq_distribution_6;
                delete _eq_distribution_7;
                delete _eq_distribution_8;
            }

            void do_preprocessing()
            {
                CONTEXT("When performing LABSWE preprocessing.");
                ///Compute bottom slopes in x and y direction
                for(unsigned long i(0); i < _grid_height; ++i)
                {
                    for(unsigned long j(0); j < _grid_width; ++j)
                    {
                        ///BASIC scheme implies: only direct predecessors are involved
                        if(i > 0 && j > 0)
                        {
                            (*_d_bottom_x)(i,j) = ((*_bottom)(i,j) - (*_bottom)(i,j - 1)) / _delta_x;
                            (*_d_bottom_y)(i,j) = ((*_bottom)(i,j) - (*_bottom)(i - 1,j)) / _delta_y;
                        }
                        else if(i == 0 && j == 0)
                        {
                            (*_d_bottom_x)(i,j) = 0.;
                            (*_d_bottom_y)(i,j) = 0.;
                        }
                        else if(i == 0)
                        {
                            (*_d_bottom_x)(i,j) = ((*_bottom)(i,j) - (*_bottom)(i,j - 1)) / _delta_x;
                            (*_d_bottom_y)(i,j) = 0;
                        }
                        else
                        {
                            (*_d_bottom_x)(i,j) = 0;
                            (*_d_bottom_y)(i,j) = ((*_bottom)(i,j) - (*_bottom)(i - 1,j)) / _delta_y;
                        }
                    }
                }
            }

            /** Capsule for the solution.
             *
             **/
            void solve()
            {
                _extract();
            };

        };
}
#endif
