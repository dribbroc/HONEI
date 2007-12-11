/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**
 * \file
 *
 * Implementation of a semi-implicit solver for the 2D Shallow Water Equations.
 * Solver is arbitrary in precision. In contrast to the RelaxSolver, which can
 * be classified as a pipelined architecture, ImplicitSolver mainly consists of
 * the assembly and solution of a linear system. Due to this, mixed precision
 * can be obtained in an intuitive way, by using an iterative refinement solver
 * for the system. Therefore, you can specify the precision of the "inner"
 * workhorse by WorkPrec_ (system assembly).
 *
 * Since we have a finite differences - discetization here, we have a slightly
 * different point of view than with RelaxSolver. Either the stepsizes in the
 * two dimensions are determined by given input fields and the dimension sizes,
 * or the other way around. For the user it will be the most popular usecase
 * to simply have some inputfields and give the dimension sizes of OMEGA,
 * assuming (0,0) to be the lower left anchor of the rectangular region.
 *
 * \ingroup grplibswe
 **/

#ifndef LIBSWE_GUARD_IMPLICIT_SOLVER_HH
#define LIBSWE_GUARD_IMPLICIT_SOLVER_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libutil/tags.hh>
#include <libmath/methods.hh>
#include <libswe/boundary_types.hh>
#include <libswe/scenario.hh>
#include <iostream>
#include <libmath/interpolation.hh>
#include <libmath/conjugate_gradients.hh>
#include <libmath/jacobi.hh>

using namespace std;
using namespace methods;
using namespace swe_solvers;
using namespace boundaries;

namespace honei {
    template<typename Tag_, typename ResPrec_, typename SolverType_, typename BoundaryType_>
    class ImplicitSolver
    {
    };

    template<typename Tag_, typename ResPrec_>
    class ImplicitSolver<Tag_, ResPrec_, CG, REFLECT>
    {
        ///Private members:
        private:
            ///Our scenario.
            Scenario<ResPrec_, IMPLICIT, REFLECT>* scenario;
            unsigned long _n;

            ///Stepsize in x direction.
            ResPrec_ _delta_x;
            ///Stepsize in y direction.
            ResPrec_ _delta_y;
            ///Size of timestep.
            ResPrec_ _delta_t;

            ///Dimensions of OMEGA:
            ResPrec_ _d_width;
            ResPrec_ _d_height;
            unsigned long _grid_width;
            unsigned long _grid_height;

            ///The input- and to-be-updated - data.
            DenseMatrix<ResPrec_> * _bottom;
            DenseMatrix<ResPrec_> * _height;
            DenseMatrix<ResPrec_> * _x_veloc;
            DenseMatrix<ResPrec_> * _y_veloc;

            ///The data to work on.
            BandedMatrix<ResPrec_> * _system_matrix;
            DenseVector<ResPrec_> * _right_hand_side;
            DenseVector<ResPrec_> * _u_temp;
            DenseVector<ResPrec_> * _v_temp;

            ///The boundary maps of the scalarfields:
            DenseMatrix<ResPrec_>* _height_bound;
            DenseMatrix<ResPrec_>* _bottom_bound;
            DenseMatrix<ResPrec_>* _x_veloc_bound;
            DenseMatrix<ResPrec_>* _y_veloc_bound;


            ///The current timestep.
            unsigned int solve_time;

            /**
             * System assembly: A.
             **/
            template<typename WorkPrec_>
            void _assemble_matrix()
            {
                DenseVector<ResPrec_> uu(_system_matrix->size(), ResPrec_(0));
                DenseVector<ResPrec_> du(_system_matrix->size(), ResPrec_(0));
                DenseVector<ResPrec_> dd(_system_matrix->size(), ResPrec_(0));
                DenseVector<ResPrec_> dl(_system_matrix->size(), ResPrec_(0));
                DenseVector<ResPrec_> ll(_system_matrix->size(), ResPrec_(0));

                ResPrec_ alpha(_delta_t*_delta_t*ResPrec_(9.81));
                ///Ignore ghost cells, assemble:
                unsigned long actual_row(1);
                unsigned long actual_column(1);
                unsigned long i(_grid_width + 3);

                while(i < (_system_matrix->size()) - (_grid_width + 2))
                {
                    ///Assemble dd:
                    dd[i] = ResPrec_(1) + ResPrec_(2)*alpha*((*_height_bound)[actual_row][actual_column])*(ResPrec_(1)/(_delta_y*_delta_y) + ResPrec_(1)/(_delta_x*_delta_x));
                    ///Assemble ll:
                    ll[i] = alpha*(((*_bottom_bound)[actual_row-1][actual_column] - (*_bottom_bound)[actual_row + 1][actual_column])/(4*_delta_y*_delta_y) - (*_height_bound)[actual_row][actual_column]/(_delta_y * _delta_y));
                    ///Assemble uu:
                    uu[i] = alpha*(((*_bottom_bound)[actual_row + 1][actual_column] - (*_bottom_bound)[actual_row - 1][actual_column])/(4*_delta_y*_delta_y) - (*_height_bound)[actual_row][actual_column]/(_delta_y * _delta_y));
                    ///Assemble dl:
                    dl[i] = alpha*(((*_bottom_bound)[actual_row][actual_column - 1] - (*_bottom_bound)[actual_row][actual_column + 1])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/(_delta_x * _delta_x));
                    ///Assemble du:
                    du[i] = alpha*(((*_bottom_bound)[actual_row][actual_column + 1] - (*_bottom_bound)[actual_row][actual_column - 1])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/(_delta_x * _delta_x));
                    ///Iterate:
                    if((actual_column) == _grid_width)
                    {
                        ++actual_row;
                        actual_column = 1;
                        i += 3;
                    }
                    else
                    {
                        ++actual_column;
                        ++i;
                    }
                }
                ///Correct boundaries:
                /*unsigned long a(0);
                unsigned long column_count(0);
                while(a < (_grid_width + 2) * (_grid_height + 2))
                {
                    cout<<a<<endl;
                    if(a < _grid_width + 2)
                    {
                        dd[a] = dd[a + _grid_width + 2];
                        du[a] = du[a + _grid_width + 2];
                        dl[a] = dl[a + _grid_width + 2];
                        uu[a] = uu[a + _grid_width + 2];
                        ll[a] = ll[a + _grid_width + 2];
                        cout<< "Accessing " << a + _grid_width + 2 << endl;
                    }

                    if(column_count == _grid_width + 1)
                    {
                        dd[a] = dd[a - 1];
                        du[a] = du[a - 1];
                        dl[a] = dl[a - 1];
                        uu[a] = uu[a - 1];
                        ll[a] = ll[a - 1];
                        cout<< "Accessing " << a -1 << endl;

                        column_count = 0;
                    }
                    else if(column_count == 0)
                    {
                        dd[a] = dd[a + 1];
                        du[a] = du[a + 1];
                        dl[a] = dl[a + 1];
                        uu[a] = uu[a + 1];
                        ll[a] = ll[a + 1];
                        cout<< "Accessing " << a + 1 << endl;

                        ++column_count;
                    }
                    else
                    {
                        ++column_count;
                    }

                    if(a > (_grid_height + 1) * (_grid_width + 2))
                    {
                        dd[a] = dd[a - _grid_width - 2];
                        du[a] = du[a - _grid_width - 2];
                        dl[a] = dl[a - _grid_width - 2];
                        uu[a] = uu[a - _grid_width - 2];
                        ll[a] = ll[a - _grid_width - 2];
                        cout<< "Accessing " << a - _grid_width - 2 << endl;

                    }

                    ++a;
                }*/
                ///Insert bands:
                _system_matrix->insert_band(0, dd);
                _system_matrix->insert_band(1, du);
                _system_matrix->insert_band(-1, dl);
                _system_matrix->insert_band(_grid_width + 2, uu);
                _system_matrix->insert_band(-_grid_width - 2, ll);

            }

            /**
             * System assembly: b.
             **/
            template<typename WorkPrec_>
            void _assemble_right_hand_side()
            {
                DenseVector<WorkPrec_> h((_grid_width +2) * (_grid_height + 2), WorkPrec_(0));

                ///Compute propagation of velocity and height fields:
                unsigned long i(1);
                unsigned long j(1);
                unsigned long current_index(_grid_width + 3);
                //\TODO: search for corruption in this loop or in interpolation:
                for(; current_index < (_grid_height + 2) * (_grid_width + 2) - (_grid_width + 2); ++current_index)
                {

                    WorkPrec_ current_x = WorkPrec_(( j - 1 ) * _delta_x);
                    WorkPrec_ current_y = WorkPrec_(( i - 1) * _delta_y);
                    h[current_index] = Interpolation<Tag_, interpolation_methods::LINEAR>::value(_delta_x, _delta_y, *_height_bound,
                            current_x - (_delta_t * (*_x_veloc_bound)[i][j]),
                            current_y - (_delta_t * (*_y_veloc_bound)[i][j]));
                    (*_u_temp)[current_index] = Interpolation<Tag_, interpolation_methods::LINEAR>::value(_delta_x, _delta_y, *_x_veloc_bound,
                            current_x - (_delta_t * (*_x_veloc_bound)[i][j]),
                            current_y - (_delta_t * (*_y_veloc_bound)[i][j]));
                    (*_v_temp)[current_index] = Interpolation<Tag_, interpolation_methods::LINEAR>::value(_delta_x, _delta_y, *_y_veloc_bound,
                            current_x - (_delta_t * (*_x_veloc_bound)[i][j]),
                            current_y - (_delta_t * (*_y_veloc_bound)[i][j]));

                    if( j == _grid_width)
                    {
                        ++i;
                        j = 1;
                        current_index += 2; //due to for
                    }
                    else
                    {
                        ++j;
                    }
                }

                ///Assemble vector:
                unsigned long current_row(1);
                unsigned long current_column(1);
                unsigned long index(_grid_width + 3);
                WorkPrec_ beta_x(_delta_t * WorkPrec_(1)/(WorkPrec_(2) * _delta_x));
                WorkPrec_ beta_y(_delta_t * WorkPrec_(1)/(WorkPrec_(2) * _delta_y));
                for(; index < (_grid_width + 2) * (_grid_height + 2) - (_grid_width + 2); ++index)
                {
                    WorkPrec_ b_diff_1((*_bottom_bound)[current_row + 1][current_column] - (*_bottom_bound)[current_row -1][current_column]);
                    WorkPrec_ b_diff_2((*_bottom_bound)[current_row][current_column + 1] - (*_bottom_bound)[current_row][current_column - 1]);
     WorkPrec_ v_x_diff, v_y_diff;

                    if(index - (_grid_width + 2) >= 0 && index + (_grid_width + 2) <= (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_x_diff = ((*_v_temp)[index + _grid_width + 2] - (*_v_temp)[index - (_grid_width + 2)]);
                    }
                    else if ( index - (_grid_width + 2) < 0 && index + (_grid_width + 2) <= (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_x_diff = ((*_v_temp)[index + _grid_width + 2] - (*_v_temp)[index]);

                    }
                    else if (index - (_grid_width + 2) >= 0 && index + (_grid_width + 2) > (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_x_diff = ((*_v_temp)[index] - (*_v_temp)[index -( _grid_width + 2)]);

                    }
                    else
                    {
                        v_x_diff = ((*_v_temp)[index] - (*_v_temp)[index]);
                    }

                    if(index - 1 >= 0 && index + 1 <= (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_y_diff = ((*_u_temp)[index + 1] - (*_u_temp)[index - 1]);
                    }
                    else if ( index - 1 < 0 && index + 1 <= (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_y_diff = ((*_u_temp)[index + 1] - (*_u_temp)[index]);

                    }
                    else if (index - 1 >= 0 && index + 1 > (_grid_width + 2) * (_grid_height + 2))
                    {
                        v_y_diff = ((*_u_temp)[index] - (*_u_temp)[index -( _grid_width + 2)]);

                    }
                    else
                    {
                        v_y_diff = ((*_u_temp)[index] - (*_u_temp)[index]);
                    }

                    //scale:
                    b_diff_1 = b_diff_1 * beta_y * (*_v_temp)[index];
                    b_diff_2 = b_diff_2 * beta_x * (*_u_temp)[index];
                    v_x_diff = v_x_diff * beta_x * (*_height_bound)[current_row][current_column];
                    v_y_diff = v_y_diff * beta_y * (*_height_bound)[current_row][current_column];
                    //accumulate:
                    (*_right_hand_side)[index] = h[index] + b_diff_1 + b_diff_2 - v_x_diff - v_y_diff;
                    //Iterate:
                    if( current_column == _grid_width)
                    {
                        ++current_row;
                        current_column = 1;
                        index += 2; //due to for
                    }
                    else
                    {
                        ++current_column;
                    }
                }

                ///Correct boundaries:
                /*unsigned long a(0);
                unsigned long column_count(0);
                while(a < (_grid_width + 2) * (_grid_height + 2))
                {
                    cout<<a<<endl;
                    if(a < _grid_width + 2)
                    {
                        (*_right_hand_side)[a] = (*_right_hand_side)[a + _grid_width + 2];
                    }

                    if(column_count == _grid_width + 1)
                    {
                        (*_right_hand_side)[a] = (*_right_hand_side)[a - 1];

                        column_count = 0;
                    }
                    else if(column_count == 0)
                    {
                        (*_right_hand_side)[a] = (*_right_hand_side)[a + 1];

                        ++column_count;
                    }
                    else
                    {
                        ++column_count;
                    }

                    if(a > (_grid_height + 1) * (_grid_width + 2))
                    {
                        (*_right_hand_side)[a] = (*_right_hand_side)[a - _grid_width - 2];

                    }

                    ++a;
                }*/


            }
            /**
             * The update of the velocity fields per timestep - also creates h.
             * */
            template<typename WorkPrec_>
            void _update(DenseVector<WorkPrec_>& w_new)
            {

                unsigned long index(_grid_width + 3);
                unsigned long actual_row(1);
                unsigned long actual_column(1);
                WorkPrec_ gamma = - WorkPrec_(9.81) * _delta_t;
                while(index < (_grid_width + 2) * (_grid_height + 2) - (_grid_width + 2))
                {
                    WorkPrec_ h_new(w_new[index] - (*_bottom_bound)[actual_row][actual_column]);

                    WorkPrec_ delta_h_1, delta_h_2;

                    if( index - 1 >= 0)
                    {
                        delta_h_1 = (h_new - (w_new[index - 1] - (*_bottom_bound)[actual_row][actual_column - 1]));
                        delta_h_2 = (h_new - (w_new[index - 1] - (*_bottom_bound)[actual_row - 1][actual_column]));
                    }
                    else
                    {
                        delta_h_1 = (h_new - (w_new[index] - (*_bottom_bound)[actual_row][actual_column - 1]));
                        delta_h_2 = (h_new - (w_new[index] - (*_bottom_bound)[actual_row - 1][actual_column]));
                    }
                    (*_height_bound)[actual_row][actual_column] = h_new;
                    (*_x_veloc_bound)[actual_row][actual_column] = gamma * delta_h_1 / _delta_x + (*_u_temp)[index];
                    (*_y_veloc_bound)[actual_row][actual_column] = gamma * delta_h_2 / _delta_y + (*_v_temp)[index];

                    //Iterate:

                    if((actual_column) == _grid_width)
                    {
                        ++actual_row;
                        actual_column = 1;
                        index += 3;
                    }
                    else
                    {
                        ++actual_column;
                        ++index;
                    }
                }
            }

        public:
            /**
             * Constructor.
             * */
            ImplicitSolver(Scenario<ResPrec_, IMPLICIT, REFLECT> & s)
            {
                scenario = &s;
            }

            /**
             * Solution capsule for one timestep.
             **/
            void solve(unsigned long iter_numbers)
            {
                //\TODO: Complete. This is for test purpose only.
                _assemble_matrix<ResPrec_>();
                _assemble_right_hand_side<ResPrec_>();

                DenseVector<ResPrec_> w_new(ConjugateGradients<Tag_, NONE>::value(*_system_matrix, *_right_hand_side, long(iter_numbers)));

                _update(w_new);

                ///Our boundary - correction:
                for(unsigned long i = 0; i < _grid_width+2; i++)
                {
                    ///Correct first row:
                    (*_height_bound)[0][i] = (*_height_bound)[1][i];
                    (*_x_veloc_bound)[0][i] = (*_x_veloc_bound)[1][i];
                    (*_y_veloc_bound)[0][i] = (*_y_veloc_bound)[1][i];
                    ///Correct last row:
                    (*_height_bound)[(_grid_height)+1][i] = (*_height_bound)[(_grid_height)][i];
                    (*_x_veloc_bound)[(_grid_height)+1][i] = (*_x_veloc_bound)[(_grid_height)][i];
                    (*_y_veloc_bound)[(_grid_height)+1][i] = (*_y_veloc_bound)[(_grid_height)][i];

                }

                for(unsigned long i = 1; i < (_grid_height)+1; i++)
                {
                    ///Correct first column:
                    (*_height_bound)[i][0] = (*_height_bound)[i][1];
                    (*_x_veloc_bound)[i][0] = (*_x_veloc_bound)[i][1];
                    (*_y_veloc_bound)[i][0] = (*_y_veloc_bound)[i][1];
                    ///Correct last column:
                    (*_height_bound)[i][(_grid_height) + 1] = (*_height_bound)[i][(_grid_height)];
                    (*_x_veloc_bound)[i][(_grid_height) + 1] = (*_x_veloc_bound)[i][(_grid_height)];
                    (*_y_veloc_bound)[i][(_grid_height) + 1] = (*_y_veloc_bound)[i][(_grid_height)];
                }
                (*_height_bound)[0][0] = (*_height_bound)[0][1];
                (*_height_bound)[0][_grid_width + 1] = (*_height_bound)[0][_grid_width];
                (*_height_bound)[_grid_height + 1][0] = (*_height_bound)[_grid_height ][0];
                (*_height_bound)[_grid_height + 1][_grid_width + 1] = (*_height_bound)[_grid_height][_grid_width];
                (*_x_veloc_bound)[0][0] = (*_x_veloc_bound)[0][1];
                (*_x_veloc_bound)[0][_grid_width + 1] = (*_x_veloc_bound)[0][_grid_width];
                (*_x_veloc_bound)[_grid_height + 1][0] = (*_x_veloc_bound)[_grid_height][0];
                (*_x_veloc_bound)[_grid_height + 1][_grid_width + 1] = (*_x_veloc_bound)[_grid_height][_grid_width];
                (*_y_veloc_bound)[0][0] = (*_y_veloc_bound)[0][1];
                (*_y_veloc_bound)[0][_grid_width + 1] = (*_y_veloc_bound)[0][_grid_width];
                (*_y_veloc_bound)[_grid_height + 1][0] = (*_y_veloc_bound)[_grid_height][0];
                (*_y_veloc_bound)[_grid_height + 1][_grid_width + 1] = (*_y_veloc_bound)[_grid_height][_grid_width];

            }

            /**
             * Solution capsule for one timestep.
             **/
            template<typename WorkPrec_>
                void solve(double conv_rad);

            /**
             * Preprocessing.
             **/
            void do_preprocessing()
            {

                ///Copy our scenario data:
                _n = scenario->n;
                this->_delta_x = scenario->delta_x;
                _delta_y = scenario->delta_y;
                _delta_t = scenario->delta_t;
                _d_width = scenario->d_width;
                _d_height = scenario->d_height;
                _grid_width = scenario->grid_width;
                _grid_height = scenario->grid_height;

                ///Make local copies of the scalarfields;
                _bottom = scenario->bottom->copy();
                _height = scenario->height->copy();
                _x_veloc = scenario->x_veloc->copy();
                _y_veloc = scenario->y_veloc->copy();
                ///Just get the address where to store boundary maps:
                _height_bound = scenario->height_bound;
                _bottom_bound = scenario->bottom_bound;
                _x_veloc_bound = scenario->x_veloc_bound;
                _y_veloc_bound = scenario->y_veloc_bound;
                ///Just get the address where to store linear system:
                _system_matrix = scenario->system_matrix;
                _right_hand_side = scenario->right_hand_side;
                _u_temp = scenario->u_temp;
                _v_temp = scenario->v_temp;
                ///Process boundary mapping:
                for(unsigned long i = 0; i < _grid_height; i++)
                {
                    for(unsigned long j = 0; j < _grid_width; j++)
                    {
                        (*(_height_bound))[i+1][j+1] = (*(_height))[i][j];
                        (*_bottom_bound)[i+1][j+1] = (*_bottom)[i][j];
                        (*_x_veloc_bound)[i+1][j+1] = (*_x_veloc)[i][j];
                        (*_y_veloc_bound)[i+1][j+1] = (*_y_veloc)[i][j];
                    }
                }
                ///Process boundary correction for all scalarfields:
                for(unsigned long i = 0; i < _grid_width+2; i++)
                {
                    ///Correct first row:
                    (*_height_bound)[0][i] = (*_height_bound)[1][i];
                    (*_bottom_bound)[0][i] = (*_bottom_bound)[1][i];
                    (*_x_veloc_bound)[0][i] = (*_x_veloc_bound)[1][i];
                    (*_y_veloc_bound)[0][i] = (*_y_veloc_bound)[1][i];
                    ///Correct last row:
                    (*_height_bound)[(_grid_height)+1][i] = (*_height_bound)[(_grid_height)][i];
                    (*_bottom_bound)[(_grid_height)+1][i] = (*_bottom_bound)[(_grid_height)][i];
                    (*_x_veloc_bound)[(_grid_height)+1][i] = (*_x_veloc_bound)[(_grid_height)][i];
                    (*_y_veloc_bound)[(_grid_height)+1][i] = (*_y_veloc_bound)[(_grid_height)][i];

                }

                for(unsigned long i = 1; i < (_grid_height)+1; i++)
                {
                    ///Correct first column:
                    (*_height_bound)[i][0] = (*_height_bound)[i][1];
                    (*_bottom_bound)[i][0] = (*_bottom_bound)[i][1];
                    (*_x_veloc_bound)[i][0] = (*_x_veloc_bound)[i][1];
                    (*_y_veloc_bound)[i][0] = (*_y_veloc_bound)[i][1];
                    ///Correct last column:
                    (*_height_bound)[i][(_grid_height) + 1] = (*_height_bound)[i][(_grid_height)];
                    (*_bottom_bound)[i][(_grid_height) + 1] = (*_bottom_bound)[i][(_grid_height)];
                    (*_x_veloc_bound)[i][(_grid_height) + 1] = (*_x_veloc_bound)[i][(_grid_height)];
                    (*_y_veloc_bound)[i][(_grid_height) + 1] = (*_y_veloc_bound)[i][(_grid_height)];
                }
                //The rest:
                (*_height_bound)[0][0] = (*_height_bound)[0][1];
                (*_height_bound)[0][_grid_width + 1] = (*_height_bound)[0][_grid_width];
                (*_height_bound)[_grid_height + 1][0] = (*_height_bound)[_grid_height ][0];
                (*_height_bound)[_grid_height + 1][_grid_width + 1] = (*_height_bound)[_grid_height][_grid_width];
                (*_bottom_bound)[0][0] = (*_bottom_bound)[0][1];
                (*_bottom_bound)[0][_grid_width + 1] = (*_bottom_bound)[0][_grid_width];
                (*_bottom_bound)[_grid_height + 1][0] = (*_bottom_bound)[_grid_height][0];
                (*_bottom_bound)[_grid_height + 1][_grid_width  + 1] = (*_bottom_bound)[_grid_height][_grid_width];
                (*_x_veloc_bound)[0][0] = (*_x_veloc_bound)[0][1];
                (*_x_veloc_bound)[0][_grid_width + 1] = (*_x_veloc_bound)[0][_grid_width];
                (*_x_veloc_bound)[_grid_height + 1][0] = (*_x_veloc_bound)[_grid_height][0];
                (*_x_veloc_bound)[_grid_height + 1][_grid_width + 1] = (*_x_veloc_bound)[_grid_height][_grid_width];
                (*_y_veloc_bound)[0][0] = (*_y_veloc_bound)[0][1];
                (*_y_veloc_bound)[0][_grid_width + 1] = (*_y_veloc_bound)[0][_grid_width];
                (*_y_veloc_bound)[_grid_height + 1][0] = (*_y_veloc_bound)[_grid_height][0];
                (*_y_veloc_bound)[_grid_height + 1][_grid_width + 1] = (*_y_veloc_bound)[_grid_height][_grid_width];

            }
    };
}
#endif
