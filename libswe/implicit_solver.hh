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
                unsigned long actual_row = 1;
                unsigned long actual_column = 1;
                unsigned long i = 0;
                while(i < (_system_matrix->size()))
                {
std::cout<<actual_row<<","<<actual_column<<endl;

                    ///Assemble dd:
                    dd[i] = ResPrec_(1) + ResPrec_(2)*alpha*((*_height_bound)[actual_row][actual_column])*(ResPrec_(1)/(_delta_x*_delta_x) + ResPrec_(1)/(_delta_y*_delta_y));
                    ///Assemble ll:
                    ll[i] = alpha*(((*_bottom_bound)[actual_row-1][actual_column] - (*_bottom_bound)[actual_row + 1][actual_column])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/_delta_x);
                    ///Assemble uu:
                    uu[i] = alpha*(((*_bottom_bound)[actual_row + 1][actual_column] - (*_bottom_bound)[actual_row - 1][actual_column])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/_delta_x);
                    ///Assemble dl:
                    dl[i] = alpha*(((*_bottom_bound)[actual_row][actual_column - 1] - (*_bottom_bound)[actual_row][actual_column + 1])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/_delta_x);
                    ///Assemble du:
                    du[i] = alpha*(((*_bottom_bound)[actual_row][actual_column + 1] - (*_bottom_bound)[actual_row][actual_column - 1])/(4*_delta_x*_delta_x) - (*_height_bound)[actual_row][actual_column]/_delta_x);
                    ///Iterate:
                    ++i;
                    if((actual_column) == _grid_width)
                    {
                        if(actual_row < _grid_height)
                        {
                            ++actual_row;
                        }
                        actual_column = 1;
                    }
                    else
                    {
                        ++actual_column;
                    }
                }

                ///Insert bands:
                _system_matrix->insert_band(0, dd);
                _system_matrix->insert_band(1, du);
                _system_matrix->insert_band(-1, dl);
                _system_matrix->insert_band(_grid_width, uu);
                _system_matrix->insert_band(-_grid_width, ll);

            }

            /**
             * System assembly: b.
             **/
            template<typename WorkPrec_>
            void _assemble_right_hand_side()
            {
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
                _delta_x = scenario->delta_x;
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
                    (*_height_bound)[(_grid_height)+1][i] = (*_height_bound)[(_grid_height) - 1][i];
                    (*_bottom_bound)[(_grid_height)+1][i] = (*_bottom_bound)[(_grid_height) - 1][i];
                    (*_x_veloc_bound)[(_grid_height)+1][i] = (*_x_veloc_bound)[(_grid_height) - 1][i];
                    (*_y_veloc_bound)[(_grid_height)+1][i] = (*_y_veloc_bound)[(_grid_height) - 1][i];

                }

                for(unsigned long i = 1; i < (scenario->grid_height)+1; i++)
                {
                    ///Correct first column:
                    (*scenario->height_bound)[i][0] = (*scenario->height_bound)[i][1];
                    (*scenario->bottom_bound)[i][0] = (*scenario->bottom_bound)[i][1];
                    (*scenario->x_veloc_bound)[i][0] = (*scenario->x_veloc_bound)[i][1];
                    (*scenario->y_veloc_bound)[i][0] = (*scenario->y_veloc_bound)[i][1];
                    ///Correct last column:
                    (*scenario->height_bound)[i][(scenario->grid_height) + 1] = (*scenario->height_bound)[i][(scenario->grid_height)];
                    (*scenario->bottom_bound)[i][(scenario->grid_height) + 1] = (*scenario->bottom_bound)[i][(scenario->grid_height)];
                    (*scenario->x_veloc_bound)[i][(scenario->grid_height) + 1] = (*scenario->x_veloc_bound)[i][(scenario->grid_height)];
                    (*scenario->y_veloc_bound)[i][(scenario->grid_height) + 1] = (*scenario->y_veloc_bound)[i][(scenario->grid_height)];
                }
            }
    };
}
#endif
