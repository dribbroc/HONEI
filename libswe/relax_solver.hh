/* vim: set number sw=4 sts=4 et nofoldenable : */

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


#ifndef LIBSWE_GUARD_RELAX_SOLVER_HH
#define LIBSWE_GUARD_RELAX_SOLVER_HH 1

/**
 * \file
 *
 * DECLARATION of a Relaxation - based, Finite Volume MUSCL in space,
 * implicit/explicit Runge Kutta in time solver. Under heavy construction.
 *
 * Solver is arbitrary in precision.
 *
 * The computed solution will converge to the unrelaxed solution (eps->0), if the
 * following conditions are satisfied:
 *
 * 1. CFL: MAX((MAXc_i)*_delta_t/_delta_x,(MAXd_i)*_delta_t/_delta_y) <=1/2 ;
 *
 * 2. Subcharacteristics: (lambda_i)²/(c_i)² + (mu_i)²*(d_i)² <= 1 ;
 *    Where lambda and mu are the eigenvalues of the Jacobian dF(u)/du or dG(u)/du
 *    respectively.
 *
 * 3. _delta_t >> _eps;
 *
 * 4. 0 < _eps << 1.
 *
 * \ingroup grplibswe
 **/

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/banded_matrix.hh>
#include <cmath>
#include <libswe/limiter.hh>
#include <libla/scale.hh>
#include <libla/dot_product.hh>
#include <libla/sum.hh>
#include <libla/element_product.hh>
#include <libla/product.hh>
#include <libutil/tags.hh>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <libswe/flow_processing.hh>
#include <libswe/source_processing.hh>
#include <libswe/post_processing.hh>
#include <libswe/scenario.hh>
#include <libswe/assembly_processing.hh>
#include <libswe/correction_processing.hh>

using namespace std;
using namespace honei;

using namespace boundaries;
using namespace directions;
using namespace assembly_types;
using namespace source_types;
using namespace swe_solvers;
namespace honei {

    typedef unsigned long ulint;
    typedef unsigned int uint;
    typedef unsigned short usint;

    template<typename Tag_,
             typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_,
             typename SourceType_,
             typename BoundaryType_>
    class RelaxSolver
    {
    };

    template<typename Tag_,
             typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_,
             typename SourceType_>
    class RelaxSolver<Tag_, ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_, SourceType_, REFLECT>
    {
        ///Private members.
        private:

            ///Our scenario:
            Scenario<ResPrec_, RELAX, REFLECT> * scenario;

            ///Stepsize in x direction.
            ResPrec_ _delta_x;
            ///Stepsize in y direction.
            ResPrec_ _delta_y;
            ///Size of timestep.
            ResPrec_ _delta_t;
            ///Current timestep.
            ulint _solve_time;

            ///The input- and to-be-updated - data.
            DenseMatrix<ResPrec_> * _bottom;
            DenseMatrix<ResPrec_> * _height;
            DenseMatrix<ResPrec_> * _x_veloc;
            DenseMatrix<ResPrec_> * _y_veloc;

            ///The relaxation parameter: 0<eps<<1.
            double _eps;

            ///The squared Manning-Coefficient used by the sourceterm computation.
            ResPrec_ _manning_n_squared;

            ///The number of cells in the finite volume descretization grid.
            ulint _n;

            ///Dimension sizes for rectangular grids.
            ulint _d_width;
            ulint _d_height;

            ///Vector _c is the relaxation - Matrix C`s diagonal-vector (must be 3-dim.).
            DenseVector<ResPrec_> * _c;

            ///Vector _d is the relaxation - Matrix D`s diagonal-vector (must be 3-dim.).
            DenseVector<ResPrec_> * _d;

            ///Vectors _u, _v, _w pointers are the relaxation vectors. size is 3N, where N is the total number of grid cells.
            ///If using boundary-mapping, the size is 3N + 4(w + h + 4).
            DenseVector<ResPrec_> * _u;
            DenseVector<ResPrec_> * _v;
            DenseVector<ResPrec_> * _w;
            DenseVector<ResPrec_> * _u_temp;
            DenseVector<ResPrec_> * _v_temp;
            DenseVector<ResPrec_> * _w_temp;


            ///Vectors for the bottom slopes.
            DenseVector<ResPrec_> * _bottom_slopes_x;
            DenseVector<ResPrec_> * _bottom_slopes_y;

            /** Encapsulates the linear combination for prediction. Uses source()
              * as well as the matrix assembling procedures to create the banded matrices.
              * Is used by solve().
              *
              * \param predictedu The temporary u vector to work on.
              * \param predictedv The temporary v vector to work on.
              * \param predictedw The temporary w vector to work on.
              *
              **/
            template<typename WorkPrec_>
            void _do_prediction(DenseVector<WorkPrec_>& predictedu, DenseVector<WorkPrec_>& predictedv, DenseVector<WorkPrec_>& predictedw)
            {
                BandedMatrix<WorkPrec_> m1(_u->size());
                BandedMatrix<WorkPrec_> m2(_u->size());
                BandedMatrix<WorkPrec_> m3(_u->size());
                BandedMatrix<WorkPrec_> m4(_u->size());
#ifdef SOLVER_BENCHMARK
                timeval start, end;
                gettimeofday(&start, 0);
#endif

                AssemblyProcessing<tags::CPU, MAIN::M1M3>::value(m1, m3, predictedu, predictedv, _delta_t, _delta_x, _d_width, _d_height, *_c);
                AssemblyProcessing<tags::CPU, MAIN::M2M4>::value(m2, m4, predictedu, predictedw, _delta_t, _delta_y, _d_width, _d_height, *_d);

                BandedMatrix<WorkPrec_> m5c(m3.copy());

                BandedMatrix<WorkPrec_> m6(_u->size());
                AssemblyProcessing<tags::CPU, QUICK::M6>::value(m1, m6, *_c, _d_width, _d_height);

                BandedMatrix<WorkPrec_> m7c(m4.copy());

                BandedMatrix<WorkPrec_> m8(_u->size());
                AssemblyProcessing<tags::CPU, QUICK::M8>::value(m2, m8, *_d, _d_width, _d_height);
#ifdef SOLVER_BENCHMARK
                gettimeofday(&end, 0);
                cout << "Assembly TOE: "<< (end.tv_sec - start.tv_sec) << " " << (end.tv_usec - start.tv_usec)<< endl;
#endif

                DenseVector<WorkPrec_> temp_u_c(predictedu.copy());
                DenseVector<WorkPrec_> temp_v_c(predictedv.copy());
                DenseVector<WorkPrec_> temp_u2_c(predictedu.copy());
                DenseVector<WorkPrec_> temp_w_c(predictedw.copy());

#ifdef SOLVER_BENCHMARK
                timeval start1, start2, end1, end2;
                gettimeofday(&start1, 0);
#endif
                DenseVector<WorkPrec_> temp1 = Product<Tag_>::value(m1,temp_v_c);
                DenseVector<WorkPrec_> temp2 = Product<Tag_>::value(m2,temp_w_c);

                DenseVector<WorkPrec_> temp3 = Product<Tag_>::value(m3,temp_u_c);
                DenseVector<WorkPrec_> temp4 = Product<Tag_>::value(m4,temp_u2_c);
#ifdef SOLVER_BENCHMARK
                gettimeofday(&end1, 0);
#endif
                DenseVector<WorkPrec_> source_c(predictedu.copy());

                SourceProcessing<source_types::SIMPLE, tags::CPU>::value(source_c, *_bottom_slopes_x, *_bottom_slopes_y, _manning_n_squared);

                DenseVector<WorkPrec_> predicted_u_temp_c(predictedu.copy());

                Sum<Tag_>::value(predicted_u_temp_c, temp1);
                Sum<Tag_>::value(predicted_u_temp_c, temp2);
                Sum<Tag_>::value(predicted_u_temp_c, temp3);
                Sum<Tag_>::value(predicted_u_temp_c, temp4);
                Sum<Tag_>::value(predicted_u_temp_c, source_c);
#ifdef SOLVER_VERBOSE
                cout << "First accu solved.\n";
#endif
                DenseVector<WorkPrec_> temp_u3_c(predictedu.copy());
                DenseVector<WorkPrec_> temp_v2_c(predictedv.copy());
                DenseVector<WorkPrec_> temp_w2_c(predictedw.copy());
                DenseVector<WorkPrec_> temp_u4_c(predictedu.copy());
#ifdef SOLVER_BENCHMARK
                gettimeofday(&start2, 0);
#endif

                DenseVector<WorkPrec_> temp11 = Product<Tag_>::value(m5c,temp_v2_c);
                DenseVector<WorkPrec_> temp22 = Product<Tag_>::value(m6, temp_u3_c);
                DenseVector<WorkPrec_> temp33 = Product<Tag_>::value(m7c,temp_w2_c);
                DenseVector<WorkPrec_> temp44 = Product<Tag_>::value(m8, temp_u4_c);
#ifdef SOLVER_BENCHMARK
                gettimeofday(&end2, 0);
                cout << "Product TOE: "<< (end1.tv_sec - start1.tv_sec) + (end2.tv_sec - start2.tv_sec)<< " " << (end1.tv_usec - start1.tv_usec) + (end2.tv_usec - start2.tv_usec)<< endl;
#endif
                DenseVector<WorkPrec_> v_c(predictedv.copy());
                DenseVector<WorkPrec_> w_c(predictedw.copy());

                Sum<Tag_>::value(v_c, temp11);
                Sum<Tag_>::value(v_c, temp22);
                Sum<Tag_>::value(w_c, temp33);
                Sum<Tag_>::value(w_c, temp44);


                predictedu = predicted_u_temp_c.copy();
                predictedv = v_c.copy();
                predictedw = w_c.copy();

#ifdef SOLVER_VERBOSE
                cout << "Second accu solved.\n";
#endif
                //Correction of reflective boundaries
                for (unsigned long j = 0; j< 3*_d_width; ++j)
                {
                    predictedu [(j+6)] = predictedu [(j+6+9*(_d_width+4))];
                    predictedv [(j+6)] = predictedv [(j+6+9*(_d_width+4))];
                    predictedw [(j+6)] = predictedw [(j+6+9*(_d_width+4))];
                    predictedu [(3*(_d_width+4)+j+6)] = predictedu [(j+6+6*(_d_width+4))];
                    predictedv [(3*(_d_width+4)+j+6)] = predictedv [(j+6+6*(_d_width+4))];
                    predictedw [(3*(_d_width+4)+j+6)] = predictedw [(j+6+6*(_d_width+4))];
                }

                for (unsigned long i = 0; i< _d_height; ++i)
                {
                    for (unsigned long k =0; k <3; ++k)
                    {
                        predictedu [((i+2)*3*(_d_width+4)+k)] = predictedu [((i+2)*3*(_d_width+4)+9+k)];
                        predictedv [((i+2)*3*(_d_width+4)+k)] = predictedv [((i+2)*3*(_d_width+4)+9+k)];
                        predictedw [((i+2)*3*(_d_width+4)+k)] = predictedw [((i+2)*3*(_d_width+4)+9+k)];
                        predictedu [((i+2)*3*(_d_width+4)+3+k)] = predictedu [((i+2)*3*(_d_width+4)+6+k)];
                        predictedv [((i+2)*3*(_d_width+4)+3+k)] = predictedv [((i+2)*3*(_d_width+4)+6+k)];
                        predictedw [((i+2)*3*(_d_width+4)+3+k)] = predictedw [((i+2)*3*(_d_width+4)+6+k)];
                        predictedu [((i+3)*3*(_d_width+4)-3+k)] = predictedu [((i+3)*3*(_d_width+4)-12+k)];
                        predictedv [((i+3)*3*(_d_width+4)-3+k)] = predictedv [((i+3)*3*(_d_width+4)-12+k)];
                        predictedw [((i+3)*3*(_d_width+4)-3+k)] = predictedw [((i+3)*3*(_d_width+4)-12+k)];
                        predictedu [((i+3)*3*(_d_width+4)-6+k)] = predictedu [((i+3)*3*(_d_width+4)-9+k)];
                        predictedv [((i+3)*3*(_d_width+4)-6+k)] = predictedv [((i+3)*3*(_d_width+4)-9+k)];
                        predictedw [((i+3)*3*(_d_width+4)-6+k)] = predictedw [((i+3)*3*(_d_width+4)-9+k)];
                    }
                }

                for (unsigned long j=0; j< 3*_d_width; ++j)
                {
                    predictedu [((_d_height+3)*3*(_d_width+4)+6+j)] = predictedu [(_d_height*3*(_d_width+4)+6+j)];
                    predictedv [((_d_height+3)*3*(_d_width+4)+6+j)] = predictedv [(_d_height*3*(_d_width+4)+6+j)];
                    predictedw [((_d_height+3)*3*(_d_width+4)+6+j)] = predictedw [(_d_height*3*(_d_width+4)+6+j)];
                    predictedu [((_d_height+2)*3*(_d_width+4)+6+j)] = predictedu [((_d_height+1)*3*(_d_width+4)+6+j)];
                    predictedv [((_d_height+2)*3*(_d_width+4)+6+j)] = predictedv [((_d_height+1)*3*(_d_width+4)+6+j)];
                    predictedw [((_d_height+2)*3*(_d_width+4)+6+j)] = predictedw [((_d_height+1)*3*(_d_width+4)+6+j)];
                }

#ifdef SOLVER_VERBOSE
                std::cout << "Finished Prediction.\n";
#endif

            }

            /** Encapsulates the setup for the values utemp, vtemp, wtemp.
              * Uses flow  - computations.
              *
              **/
            template<typename WorkPrec_>
            void _do_setup_stage1(DenseVector<WorkPrec_>& su, DenseVector<WorkPrec_>& sv, DenseVector<WorkPrec_>& sw)
            {
                WorkPrec_ prefac(0);
                if(_eps != _delta_t)
                {
                    prefac = WorkPrec_(1)/(_eps - _delta_t);
                }
                else
                {
                    cout << "prefac is invalid!\n";
                }
                DenseVector<WorkPrec_> vc(_v->copy());
                DenseVector<WorkPrec_> u1_c(_u->copy());
                DenseVector<WorkPrec_> u2_c(_u->copy());
                FlowProcessing<directions::X, tags::CPU>::value(u1_c);

                Scale<Tag_>::value(vc, _eps);
                Scale<Tag_>::value(u1_c, -_delta_t);
                Sum<Tag_>::value(vc, u1_c);
                DenseVector<WorkPrec_> tempsum(vc.copy());

                Scale<Tag_>::value(tempsum, prefac);
                DenseVector<WorkPrec_> wc(_w->copy());
                FlowProcessing<directions::Y, tags::CPU>::value(u2_c);

                Scale<Tag_>::value(wc, _eps);
                Scale<Tag_>::value(u2_c, -_delta_t);
                Sum<Tag_>::value(wc, u2_c);
                DenseVector<WorkPrec_> tempsum2(wc.copy());

                Scale<Tag_>::value(tempsum2, prefac);
#ifdef SOLVER_VERBOSE
                cout << "Temp relax vectors after building:\n";
                cout << stringify(*_u_temp) << endl;
                cout << stringify(*_v_temp) << endl;
                cout << stringify(*_w_temp) << endl;
                std::cout << "Finished Setup 1.\n";
#endif
                sv = tempsum.copy();
                sw = tempsum2.copy();

            }

            /** Encapsulates the setup for the values utemp, vtemp, wtemp.
              * Uses flow - computations.
              *
              * \param predictedu The vector u to work on.
              * \param predictedv The vector v to work on.
              * \param predictedw The vector w to work on.
              *
              **/
            template<typename WorkPrec_>
            void _do_setup_stage2(DenseVector<WorkPrec_>& predictedu,
                                  DenseVector<WorkPrec_>& predictedv,
                                  DenseVector<WorkPrec_>& predictedw)
            {
                ///Apply flow to newest u:
                DenseVector<WorkPrec_> f_c(predictedu.copy());
                FlowProcessing<directions::X, tags::CPU>::value(f_c);

                ///Compute linear combinations und accumulate:
                DenseVector<WorkPrec_> v_result_c(predictedv.copy());
                Scale<Tag_>::value(v_result_c, _eps);
                Scale<Tag_>::value(f_c, _delta_t);
                Sum<Tag_>::value(f_c, v_result_c);
                DenseVector<WorkPrec_> innersum1(v_result_c.copy());

                ///Apply flow to old u:

                DenseVector<WorkPrec_> flow_c(_u_temp->copy());
                FlowProcessing<directions::X, tags::CPU>::value(flow_c);

                DenseVector<WorkPrec_> v_temp_result_c(_v_temp->copy());

                Scale<Tag_>::value(v_temp_result_c, WorkPrec_(-2)*_delta_t);
                Scale<Tag_>::value(flow_c, WorkPrec_(2)*_delta_t);
                Sum<Tag_>::value(v_temp_result_c, flow_c);
                DenseVector<WorkPrec_> innersum2(v_temp_result_c.copy());

                Sum<Tag_>::value(innersum1, innersum2);

                ///Scale sum:
                Scale<Tag_>::value(innersum1, WorkPrec_(1)/(_eps+_delta_t));

                ///Repeat for w:
                DenseVector<WorkPrec_> flow2_c(predictedu.copy());
                FlowProcessing<directions::Y, tags::CPU>::value(flow2_c);

                DenseVector<WorkPrec_> w_result_c(predictedw.copy());

                Scale<Tag_>::value(w_result_c, _eps);
                Scale<Tag_>::value(flow2_c, _delta_t);
                Sum<Tag_>::value(flow2_c, w_result_c);
                DenseVector<WorkPrec_> innersum11(w_result_c.copy());

                DenseVector<WorkPrec_> flow3_c(_u_temp->copy());

                FlowProcessing<directions::Y, tags::CPU>::value(flow3_c);

                DenseVector<WorkPrec_> w_temp_result_c(_w_temp->copy());

                Scale<Tag_>::value(w_temp_result_c, WorkPrec_(-2)*_delta_t);
                Scale<Tag_>::value(flow3_c, WorkPrec_(2)*_delta_t);
                Sum<Tag_>::value(w_temp_result_c, flow3_c);
                DenseVector<WorkPrec_>innersum22(w_temp_result_c.copy());

                Sum<Tag_>::value(innersum11, innersum22);
                Scale<Tag_>::value(innersum11, WorkPrec_(1)/(_eps + _delta_t));

                predictedv = innersum1.copy();
                predictedw = innersum11.copy();
#ifdef SOLVER_VERBOSE
                std::cout << "Finished Setup 2.\n";
#endif
            }

        public:
            /**
              * Returns the current renderable heigth matrix
              *
              **/
            DenseMatrix<ResPrec_> &getHeight()
            {
                return *_height;
            }

            /**
              * Returns the renderable bottom.
              *
              **/
            DenseMatrix<ResPrec_> &getBottom()
            {
                return *_bottom;
            }

            /**
              * Performs the preprocessing.
              *
              * Implementation of the preprocessing stage of the RelaxSolver.
              *
              * At first, the input scalarfields have to be mapped onto a larger
              * DenseMatrix in order to apply boundary conditions.
              *
              * Secondly, the h,u1 und u2 - values have to be written to the
              * locations in the relaxation vectors.
              *
              * The preprocessing stage`s third task is to compute the bottom slopes.
              **/
            void do_preprocessing()
            {
                /// Setting up initial conditions:
                /// v_0(x,y) = F(u_0(x,y)) using flowX(),
                /// w_0(x,y) = G(u_0(x,y)) using flowY().
                /// where u = transpose((h q1 q2)), v = transpose((v1 v2 v3)) and w = transpose((w1 w2 w3))
                /// and q1 = h*u1, q2 = h*u2.
                /// Then apply boundary conditions.

                ///Provide maps.
                DenseMatrix<ResPrec_> hbound((this->_d_height)+4,  (this->_d_width)+4, 0);
                DenseMatrix<ResPrec_> u1bound((this->_d_height)+4, (this->_d_width)+4, 0);
                DenseMatrix<ResPrec_> u2bound((this->_d_height)+4, (this->_d_width)+4, 0);
                DenseMatrix<ResPrec_> bbound((this->_d_height)+4,  (this->_d_width)+4, 0);
#ifdef  SOLVER_VERBOSE
                std::cout << "Preproc: Maps provided.\n";
#endif
                //Boundary Map
                //0 0 1 1 1 1 0 0
                //0 0 2 2 2 2 0 0
                //3 4 x x x x 5 6
                //3 4 x x x x 5 6
                //3 4 x x x x 5 6
                //3 4 x x x x 5 6
                //0 0 7 7 7 7 0 0
                //0 0 8 8 8 8 0 0
                //
                for (unsigned long j = 0; j< _d_width; ++j)
                {
                    //Setting boundary 1 to second matrix row
                    hbound[0][j+2] = (*_height)[1][j];
                    bbound[0][j+2] = (*_bottom)[1][j];
                    u1bound[0][j+2] = (*_x_veloc)[1][j];
                    u2bound[0][j+2] = (*_y_veloc)[1][j];
                    //Setting boundary 2 to first matrix row
                    hbound[1][j+2] = (*_height)[0][j];
                    bbound[1][j+2] = (*_bottom)[0][j];
                    u1bound[1][j+2] = (*_x_veloc)[0][j];
                    u2bound[1][j+2] = (*_y_veloc)[0][j];
                }

                for (unsigned long i = 0; i< _d_height; ++i)
                {
                    //Setting boundary 3 to second matrix column
                    hbound[i+2][0] = (*_height)[i][1];
                    bbound[i+2][0] = (*_bottom)[i][1];
                    u1bound[i+2][0] = (*_x_veloc)[i][1];
                    u2bound[i+2][0] = (*_y_veloc)[i][1];
                    //Setting boundary 4 to first matrix column
                    hbound[i+2][1] = (*_height)[i][0];
                    bbound[i+2][1] = (*_bottom)[i][0];
                    u1bound[i+2][1] = (*_x_veloc)[i][0];
                    u2bound[i+2][1] = (*_y_veloc)[i][0];
                    //Take over inner values
                    for(unsigned long j=0; j< _d_width; ++j)
                    {
                        hbound[i+2][j+2] = (*_height)[i][j];
                        bbound[i+2][j+2] = (*_bottom)[i][j];
                        u1bound[i+2][j+2] = (*_x_veloc)[i][j];
                        u2bound[i+2][j+2] = (*_y_veloc)[i][j];
                    }
                    //Setting boundary 5 to rightmost matrix column
                    hbound[i+2][_d_width+2] = (*_height)[i][_d_width-1];
                    bbound[i+2][_d_width+2] = (*_bottom)[i][_d_width-1];
                    u1bound[i+2][_d_width+2] = (*_x_veloc)[i][_d_width-1];
                    u2bound[i+2][_d_width+2] = (*_y_veloc)[i][_d_width-1];
                    //Setting boundary 6 to rightmost-1 matrix column
                    hbound[i+2][_d_width+3] = (*_height)[i][_d_width-2];
                    bbound[i+2][_d_width+3] = (*_bottom)[i][_d_width-2];
                    u1bound[i+2][_d_width+3] = (*_x_veloc)[i][_d_width-2];
                    u2bound[i+2][_d_width+3] = (*_y_veloc)[i][_d_width-2];
                }

                for(unsigned long j=0; j< _d_width; ++j)
                {
                    //Setting boundary 7 to last matrix row
                    hbound[_d_height+2][j+2] = (*_height)[_d_height-1][j];
                    bbound[_d_height+2][j+2] = (*_bottom)[_d_height-1][j];
                    u1bound[_d_height+2][j+2] = (*_x_veloc)[_d_height-1][j];
                    u2bound[_d_height+2][j+2] = (*_y_veloc)[_d_height-1][j];
                    //Setting boundary 8 to last-1 matrix row
                    hbound[_d_height+3][j+2] = (*_height)[_d_height-2][j];
                    bbound[_d_height+3][j+2] = (*_bottom)[_d_height-2][j];
                    u1bound[_d_height+3][j+2] = (*_x_veloc)[_d_height-2][j];
                    u2bound[_d_height+3][j+2] = (*_y_veloc)[_d_height-2][j];
                }
#ifdef SOLVER_VERBOSE
                std::cout << "Preproc: Mapping done.\n";
                cout << stringify(hbound) << endl;
                cout << stringify(bbound) << endl;
                cout << stringify(u1bound) << endl;
                cout << stringify(u2bound) << endl;
#endif

                ///Building up the relaxation - vectors by concatenating the maps` rows.
                ///We need to compute u first in order to be able to compute the initial flows. After this, by using
                ///forward iterators, the v and w vectors can be set up.
                typename DenseVector<ResPrec_>::ElementIterator k(_u->begin_elements());
                for (ulint i= 0; i!= hbound.rows(); ++i)
                {


                    DenseVectorRange<ResPrec_> actual_row = hbound[i];
                    for(typename DenseVector<ResPrec_>::ElementIterator j(actual_row.begin_elements()),
                            j_END(actual_row.end_elements());
                            //k((*_u).begin_elements());
                            j!= j_END; ++j)
                    {
                        (*_u)[k.index()] = hbound[i][j.index()];
                        (*_u)[(k.index())+1] = u1bound[i][j.index()] * hbound[i][j.index()];
                        (*_u)[(k.index())+2] = u2bound[i][j.index()] * hbound[i][j.index()];
                        ++k; ++k; ++k;
                    }
                }
#ifdef SOLVER_VERBOSE
                cout << "u^T after building:\n";
                cout << stringify(*_u) << endl;
#endif
                DenseVector<ResPrec_> uFlow = _u->copy();
                FlowProcessing<directions::X, tags::CPU>::value(uFlow);
                (*_v) = uFlow;
#ifdef SOLVER_VERBOSE
                cout << "v^T after building:\n";
                cout << stringify(*_v) << endl;
#endif

                DenseVector<ResPrec_> u2Flow = _u->copy();
                FlowProcessing<directions::Y, tags::CPU>::value(u2Flow);
                (*_w) = u2Flow;
#ifdef SOLVER_VERBOSE
                cout << "w^T after building:\n";
                cout << stringify(*_w) << endl;
#endif
                ///Now, that the relaxation vectors have been provided, the only thing left to do is to
                ///compute the bottom slopes.
                typename DenseVector<ResPrec_>::ElementIterator l(_bottom_slopes_x->begin_elements());
                typename DenseVector<ResPrec_>::ElementIterator k4(_bottom_slopes_y->begin_elements());
                for (ulint i = 0; i!= bbound.rows(); ++i)
                {
                    DenseVectorRange<ResPrec_> actual_row = bbound[i];
                    for(typename DenseVector<ResPrec_>::ConstElementIterator j(actual_row.begin_elements()),
                            j_END(actual_row.end_elements());
                            //k((*_bottom_slopes_x).begin_elements()),
                            //l((*_bottom_slopes_y).begin_elements());
                            j!= j_END; ++j)
                    {
                        if(i>0 /*&& j.index()>0*/)
                        {
                            (*_bottom_slopes_y)[k4.index()] = (bbound[i][j.index()] - bbound[i-1][j.index()]) /this->_delta_y;
                            //(*_bottom_slopes_x)[l.index()] = (bbound[i][j.index()] - bbound[i][(j.index())-1]) /this->_delta_x;
                        }
                        else
                        {
                            //(*_bottom_slopes_x)[k4.index()] = -100000;
                            (*_bottom_slopes_y)[k4.index()] = 0;
                        }
                        if(j.index()>0)
                        {
                            //(*_bottom_slopes_y)[k4.index()] = (bbound[i][j.index()] - bbound[i-1][j.index()]) /this->_delta_y;
                            (*_bottom_slopes_x)[l.index()] = (bbound[i][j.index()] - bbound[i][(j.index())-1]) /this->_delta_x;

                        }
                        else
                        {
                            (*_bottom_slopes_x)[l.index()] = 0;
                            //(*_bottom_slopes_y)[l.index()] = -100000;

                        }


                        ++k4;
                        ++l;
                    }
                }
#ifdef SOLVER_VERBOSE
                cout << "Slopes after building:\n";
                cout << stringify(*_bottom_slopes_x) << endl;
                cout << stringify(*_bottom_slopes_y) << endl;
                std::cout << "Finished preprocessing.\n";
#endif

            }

            ///Constructors
            /**
             * First simple public constructor for tests.
             *
             * \param heigth The input heigth -field.
             * \param bottom The input bottom -field.
             * \param u1 The x -velocity field.
             * \param u2 The y -velocity field.
             * \param u The relaxation vector u.
             * \param v The relaxation vector v.
             * \param w The relaxation vector w.
             * \param dwidth The width of the FV - discretization grid.
             * \param dheigth The heigth of the FV - discretization grid.
             * \param deltax The x - stepsize.
             * \param deltay The y - stepsize.
             * \param deltat The time - stepsize.
             * \param eps The relaxation parameter.
             * \param bottomx The vector for the bottom slopes in x direction.
             * \param bottomy The vector for the bottom slopes in y direction.
             * \param c The vector c.
             * \param d The vector d.
             *
             **/
            RelaxSolver(DenseMatrix<ResPrec_> *height,
                    DenseMatrix<ResPrec_> *bottom,
                    DenseMatrix<ResPrec_> *u1,
                    DenseMatrix<ResPrec_> *u2,
                    DenseVector<ResPrec_> *u,
                    DenseVector<ResPrec_> *v,
                    DenseVector<ResPrec_> *w,
                    ulint dwidth,
                    ulint dheight,
                    ResPrec_ deltax,
                    ResPrec_ deltay,
                    ResPrec_ deltat,
                    double eps,
                    DenseVector<ResPrec_> * bottomx,
                    DenseVector<ResPrec_> * bottomy,
                    DenseVector<ResPrec_> * c,
                    DenseVector<ResPrec_> * d,
                    ResPrec_ manning_n)
            {
                this->_height = height;
                this->_bottom = bottom;
                this->_x_veloc = u1;
                this->_y_veloc = u2;
                this->_u = u;
                this->_v = v;
                this->_w = w;
                this->_d_width = dwidth;
                this->_d_height = dheight;
                this->_delta_x = deltax;
                this->_delta_y = deltay;
                this->_delta_t = deltat;
                this->_eps = eps;
                this->_solve_time = 0;

                this->_n = _d_width * _d_height;

                this->_bottom_slopes_x = bottomx;
                this->_bottom_slopes_y = bottomy;

                this->_c = c;
                this->_d = d;
                this->_manning_n_squared = manning_n * manning_n;
                /*_u_temp = new DenseVector<ResPrec_> (_u->copy());
                  _v_temp = new DenseVector<ResPrec_> (_v->copy());
                  _w_temp = new DenseVector<ResPrec_> (_w->copy());
                  */
            }
            /**
             * Second simple public constructor for tests.
             *
             * \param scenario Our scenario.
             *
             **/
            RelaxSolver(Scenario<ResPrec_, swe_solvers::RELAX, REFLECT> & scenario )
            {
                this->scenario = &scenario;
                this->_height = scenario.height;
                this->_bottom = scenario.bottom;
                this->_x_veloc = scenario.x_veloc;
                this->_y_veloc = scenario.y_veloc;
                this->_u = scenario.u;
                this->_v = scenario.v;
                this->_w = scenario.w;
                this->_d_width = scenario.d_width;
                this->_d_height = scenario.d_height;
                this->_delta_x = scenario.delta_x;
                this->_delta_y = scenario.delta_y;
                this->_delta_t = scenario.delta_t;
                this->_eps = scenario.eps;
                this->_solve_time = 0;
                this->_n = _d_width * _d_height;

                this->_bottom_slopes_x = scenario.bottom_slopes_x;
                this->_bottom_slopes_y = scenario.bottom_slopes_y;

                this->_c = scenario.c;
                this->_d = scenario.d;
                this->_manning_n_squared = scenario.manning_n * scenario.manning_n;
            }

            /** Encapsulates computation in one timestep. In the driver-application, one
             * can simply write a loop in which solve is called at first and then the
             * renderable matrices are read out.
             *
             **/
            void solve()
            {
                DenseVector<ResPrec_> u_temp(_u->copy());//(_u->size());
                DenseVector<ResPrec_> v_temp(_v->copy());//(_u->size());
                DenseVector<ResPrec_> w_temp(_w->copy());//(_u->size());

                _u_temp = &u_temp;
                _v_temp = &v_temp;
                _w_temp = &w_temp;
                _do_setup_stage1<InitPrec1_>(*_u_temp, *_v_temp, *_w_temp );

                DenseVector<PredictionPrec1_> predictedu_c(_u_temp->copy());
                DenseVector<PredictionPrec1_> predictedv_c(_v_temp->copy());
                DenseVector<PredictionPrec1_> predictedw_c(_w_temp->copy());

                _do_prediction<PredictionPrec1_>(predictedu_c, predictedv_c, predictedw_c);
                _do_setup_stage2<InitPrec2_>(predictedu_c, predictedv_c, predictedw_c);

                _do_prediction<PredictionPrec2_>(predictedu_c, predictedv_c, predictedw_c);
                CorrectionProcessing<REFLECT, tags::CPU>::value(predictedu_c, predictedv_c, predictedw_c, *_u, *_v, *_w, _d_width, _d_height, *_height);
                ++_solve_time;

#ifdef SOLVER_VERBOSE
                cout << "Corrected u, finished solution, timestep:" << stringify(_solve_time) << endl;
                cout << stringify(*_u)<<endl;
#endif

#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(1);
#endif

            }
    };

}
#endif
