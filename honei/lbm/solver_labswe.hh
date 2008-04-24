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

#include <honei/liblbm/tags.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
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

                unsigned long _grid_width, _grid_height;

                DenseVector<ResPrec_>* _equilibrium_distribution;
                DenseVector<ResPrec_>* _source_x;
                DenseVector<ResPrec_>* _source_y;

                /** Global constants.
                 *
                 **/
                ResPrec_ _n_alpha, _e, _gravity, _pi;
                DenseVector<ResPrec_>* _distribution_vector_x;
                DenseVector<ResPrec_>* _distribution_vector_y;


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
                _e = _delta_x / _delta_t;
                _distribution_vector_x = new DenseVector<ResPrec_>(9ul);
                _distribution_vector_y = new DenseVector<ResPrec_>(9ul);
                _d_bottom_x = new DenseMatrix<ResPrec_>(gx, gy);
                _d_bottom_y = new DenseMatrix<ResPrec_>(gx, gy);
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
                delete _distribution_vector_x;
                delete _distribution_vector_y;
                delete _d_bottom_x;
                delete _d_bottom_y;
            }

            void do_preprocessing()
            {
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
        };
}
#endif
