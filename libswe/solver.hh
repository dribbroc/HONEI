/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Volker Jung <volker.m.jung@t-online.de>
 * Copyright (c) 2007 Joachim Messer <joachim.messer@t-online.de>
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


#ifndef LIBSWE_GUARD_SOLVER_HH
#define LIBSWE_GUARD_SOLVER_HH 1

/**
 * \file
 *
 * DECLARATION of a Relaxation - based, Finite Volume MUSCL in space,
 * implicit/explicit Runge Kutta in time solver. Under heavy construction.
 *
 * Solver is arbitrary in precision.
 *
 * \ingroup grplibswe
 **/
 
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/banded_matrix.hh>
#include <cmath>
#include <libswe/limiter.hh>
#include <libla/vector_scaled_sum.hh>
#include <libla/scalar_vector_product.hh>
#include <libla/scalar_product.hh>
#include <libla/vector_sum.hh>
#include <libla/vector_elementwise_product.hh>
#include <libla/matrix_vector_product.hh>
#include <libutil/tags.hh>
#include <iostream>
#include <fstream>

using namespace std;

namespace pg512 {

    typedef unsigned long ulint;
    typedef unsigned int uint;
    typedef unsigned short usint;

    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_>
    class RelaxSolver
    {
        ///Private members.
        private:
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
	    static const double _manning_n_squared = 0.000625;

            ///The number of cells in the finite volume descretization grid.
            ulint _n;
            
            ///Dimension sizes for rectangular grids.
            ulint _d_width;
            ulint _d_height;

            ///Vector _c is the relaxation - Matrix C`s diagonal-vector (must be 3-dim.).
            DenseVector<ResPrec_> * _c;

            ///Vector _d is the relaxation - Matrix D`s diagonal-vector (must be 3-dim.).
            DenseVector<ResPrec_> * _d;

            /** Vectors contain the boundary - scalars => They have to be (2a+2b) - dimensional,
              * if a is the number of cell-steps in x-direction and b represents the same for
              * y- direction (For square - shaped grids it is 4n dimensional, where n = squareroot(N), and N is the total number of cells)
              * This might be useful for very complicated simulation settings. In simpler cases, one should use the
              * below options.
              **/
            DenseVector<ResPrec_> * _bottom_bound;
            DenseVector<ResPrec_> * _height_bound;
            DenseVector<ResPrec_> * _xveloc_bound;
            DenseVector<ResPrec_> * _yveloc_bound;

            /**
              * Boundaries: Scalars for rectangular grids. Simple version of the above vector-
              * based solution.
              *
              **/
            ResPrec_ _north_bound;
            ResPrec_ _south_bound;
            ResPrec_ _west_bound;
            ResPrec_ _east_bound;

            ///Flags for boundary usage.
            bool _simple_bound;
            bool _usage_reflect;
            bool _usage_constant;
            bool _usage_cyclic;
            bool _usage_transmissive;
        
            ///Vectors _u, _v, _w pointers are the relaxation vectors. size is 3N, where N is the total number of grid cells.
            ///If using boundary-mapping, the size is 3N + 4(w + h + 4).
            DenseVector<ResPrec_> * _u;
            DenseVector<ResPrec_> * _v;
            DenseVector<ResPrec_> * _w;
            DenseVector<ResPrec_> * _u_temp;
            DenseVector<ResPrec_> * _v_temp;
            DenseVector<ResPrec_> * _w_temp;
            //DenseVector<ResPrec_> * _u_predicted;
            //DenseVector<ResPrec_> * _v_predicted;
            //DenseVector<ResPrec_> * _w_predicted;;



            ///Vectors for the bottom slopes.
            DenseVector<ResPrec_> * _bottom_slopes_x; //size:w
            DenseVector<ResPrec_> * _bottom_slopes_y; //size:h

            /** Basic matrix assembly. Uses Limiters and theta().
              * Computes M_1, M_3.
              *
              * \param m1 Matrix m1 is the first Matrix to be assembled.
              * \param m3 Matrix m3 is the second Matrix to be assembled.
              *
              **/
           
            template<typename WorkPrec_>
            void _assemble_matrix1(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& m3, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v);

            /** Basic matrix assembly. Uses Limiters and theta().
              * Computes M_2, M_4.
              *
              * \param m1 Matrix m2 is the first Matrix to be assembled.
              * \param m3 Matrix m4 is the second Matrix to be assembled.
              *
              **/
            template<typename WorkPrec_>
            void _assemble_matrix2(BandedMatrix<WorkPrec_>& m2, BandedMatrix<WorkPrec_>& m4, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* w);


            /** Simplified matrix assembly.
              * Computes M_5.
              *
              **/
            template<typename WorkPrec_>
            /*BandedMatrix<WorkPrec_>*/ void _quick_assemble_matrix1(BandedMatrix<WorkPrec_>& m3, BandedMatrix<WorkPrec_>& result)
	    {
                //result = *(m3.copy());
		//return ;
	    }

            /** Simplified matrix assembly.
              * Computes M_6.
              *
              **/
            template<typename WorkPrec_>
            /*BandedMatrix<WorkPrec_>*/void _quick_assemble_matrix2(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& result);

            /** Simplified matrix assembly.
              * Computes M_7.
              *
              **/
            template<typename WorkPrec_>
            /*BandedMatrix<WorkPrec_>*/void _quick_assemble_matrix3(BandedMatrix<WorkPrec_>& m4, BandedMatrix<WorkPrec_>& result)
	    {
                //result = *(m4.copy());
		//return ;
	    }

            /** Simplified matrix assembly.
              * Computes M_8.
              *
              **/
            template<typename WorkPrec_>
            /*BandedMatrix<WorkPrec_>*/void _quick_assemble_matrix4(BandedMatrix<WorkPrec_>& m2, BandedMatrix<WorkPrec_>& result);

            /** Flow computation.
              * Used by preprocessing.
              *
              * \param i Access Parameter 1.
              * \param j Access Parameter 2.
              * 
              **/
            template<typename WorkPrec_>
            DenseVector<WorkPrec_> _flow_x(uint i, uint j);

            /** Flow computation.
              * Used by preprocessing.
              *
              * \param i Access Parameter 1.
              * \param j Access Parameter 2.
              * 
              **/
            template<typename WorkPrec_>
            DenseVector<WorkPrec_> _flow_y(uint i, uint j);

            /** Flow computation.
              *
              **/
	    template<typename WorkPrec_>
            void _flow_x(DenseVector<WorkPrec_> & vector);

            /** Flow computation.
              * 
              **/
            template<typename WorkPrec_>
            void _flow_y(DenseVector<WorkPrec_> & vector);


            /** Source Term computation.
              *
              * \param i Access Parameter 1.
              * \param j Access Parameter 2.
              * 
              **/
            template<typename WorkPrec_>
            DenseVector<WorkPrec_> _source(uint i, uint j);

	    /** Source Term computation.
	     *
	     * \param vector Densevector for which the source term should be computed.
	     *
	     **/
	    template<typename WorkPrec_>
	    void _source(DenseVector<WorkPrec_>& vector);


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
            void _do_prediction(DenseVector<WorkPrec_>& predictedu, DenseVector<WorkPrec_>& predictedv, DenseVector<WorkPrec_>& predictedw);

            /** Encapsulates the setup for the values utemp, vtemp, wtemp.
              * Uses flow  - computations.
              *
              **/
            template<typename WorkPrec_>
            void _do_setup_stage1( DenseVector<WorkPrec_>& su, DenseVector<WorkPrec_>& sv, DenseVector<WorkPrec_>& sw);

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
                                  DenseVector<WorkPrec_>& predictedw);

            /** Encapsulates the correction stage.
              *  Precision is that of the result.
              * \param predu The temporary vector u.
              * \param predv The temporary vector v.
              * \param predw The temporary vector w.
             **/ 
            void _do_correction(DenseVector<ResPrec_>& predu,DenseVector<ResPrec_>& predv,DenseVector<ResPrec_>& predw);

            /** First simple postprocessing method, generating a GNUPLOT>splot
              * compatible matrix file for the height field of the current timestep.
              *
              * \param every Generate file for every xth timestep.
              **/ 
            void _do_postprocessing(int every);

            /** Basic matrix assembly. Uses Limiters. DEBUG version.
              * Computes M_1, M_3.
              *
              * \param m1 Matrix m1 is the first Matrix to be assembled.
              * \param m3 Matrix m3 is the second Matrix to be assembled.
              *
              **/
            template<typename WorkPrec_>
            void _assemble_matrix1_DEBUG(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& m3, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v);

            /** Basic matrix assembly. Uses Limiters. DEBUG version.
              * Computes M_2, M_4.
              *
              * \param m2 Matrix m1 is the first Matrix to be assembled.
              * \param m4 Matrix m3 is the second Matrix to be assembled.
              *
              **/
           
            template<typename WorkPrec_>
            void _assemble_matrix2_DEBUG(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& m3, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v);


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
              **/
            void do_preprocessing();

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
                        DenseVector<ResPrec_> * d)
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

                this->_simple_bound = true;
                this->_usage_reflect = true;
                this->_usage_constant = false;
                this->_usage_cyclic = false;
                this->_usage_transmissive = false;

                this->_n = _d_width * _d_height;

                this->_bottom_slopes_x = bottomx;
                this->_bottom_slopes_y = bottomy;

                this->_c = c;
                this->_d = d;
                _u_temp = _u->copy();//reinterpret_cast<DenseVector<WorkPrec_>*>(_u);
                _v_temp = _v->copy();//reinterpret_cast<DenseVector<WorkPrec_>*>(_v);
                _w_temp = _w->copy();//reinterpret_cast<DenseVector<WorkPrec_>*>(_w);
                //_u_predicted = new DenseVector<ResPrec_>(*(_u_temp->copy()));
                //_v_predicted = new DenseVector<ResPrec_>(*(_v_temp->copy()));
                //_w_predicted = new DenseVector<ResPrec_>(*(_w_temp->copy()));


            }
            /** Encapsulates computation in one timestep. In the driver-application, one
              * can simply write a loop in which solve is called at first and then the
              * renderable matrices are read out.
              *
              **/
            void solve();
    };

    ///MEMBER FUNCTION TEMPLATE IMPLEMENTATION

    /**
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
     *
     **/
    template<typename ResPrec_, 
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::do_preprocessing()
    {   
        /// Setting up initial conditions:
        /// v_0(x,y) = F(u_0(x,y)) using flowX(),
        /// w_0(x,y) = G(u_0(x,y)) using flowY().
        /// where u = transpose((h q1 q2)), v = transpose((v1 v2 v3)) and w = transpose((w1 w2 w3))
        /// and q1 = h*u1, q2 = h*u2.
        /// Then apply boundary conditions.
    
        ///Provide maps.
        DenseMatrix<ResPrec_> hbound((this->_d_width)+4,  (this->_d_height)+4, 0);
        DenseMatrix<ResPrec_> u1bound((this->_d_width)+4, (this->_d_height)+4, 0);
        DenseMatrix<ResPrec_> u2bound((this->_d_width)+4, (this->_d_height)+4, 0);
        DenseMatrix<ResPrec_> bbound((this->_d_width)+4,  (this->_d_height)+4, 0);
        std::cout << "Preproc: Maps provided.\n";

        ///Do the mapping by applying boundary - usage.
        //if(this->_usage_reflect && this->_simple_bound)
        //{
            ///If assuming, that all input fields are exactly of the same size, we can do all the work within
            ///one loop - pair:
            for (unsigned long i = 0; i!= hbound.rows(); ++i) 
            {
                DenseVector<ResPrec_> actual_row = hbound[i];
                for(typename DenseVector<ResPrec_>::ElementIterator j(actual_row.begin_elements()),
                                                                     j_END(actual_row.end_elements());
                                                                     j!= j_END; ++j)
                {
                    ///Check if boundary - ghost cell is going to be accessed.
                    if(i<2 || i>=(hbound.rows()-2) ||(j.index()<2 || j.index() >=(hbound.columns()-2)))
                    {
                        hbound[i][j.index()] = 5;
                        bbound[i][j.index()] = 1; 
                        u1bound[i][j.index()] = 1;
                        u2bound[i][j.index()] = 1;
                    }
                    else
                    {
                        hbound[i][j.index()] = (*_height)[i-2][(j.index())-2];
                        bbound[i][j.index()] = (*_bottom)[i-2][(j.index())-2];
                        u1bound[i][j.index()] =(*_x_veloc)[i-2][(j.index())-2];
                        u2bound[i][j.index()] =(*_y_veloc)[i-2][(j.index())-2];
                    }
                }
            }
            std::cout << "Preproc: Mapping done.\n";
            cout << stringify(hbound) << endl;
            cout << stringify(bbound) << endl;
            cout << stringify(u1bound) << endl;
            cout << stringify(u2bound) << endl;

       // }//TODO: the other cases of boundary usage

        ///Building up the relaxation - vectors by concatenating the maps` rows.
        ///We need to compute u first in order to be able to compute the initial flows. After this, by using
        ///forward iterators, the v and w vectors can be set up.
        typename DenseVector<ResPrec_>::ElementIterator k(_u->begin_elements());
        for (ulint i= 0; i!= hbound.rows(); ++i) 
        {


            DenseVector<ResPrec_> actual_row = hbound[i];
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
        cout << "u^T after building:\n";
        cout << stringify(*_u) << endl;
        /*OBSOLETE
        typename DenseVector<ResPrec_>::ElementIterator k2(_v->begin_elements());
        for (ulint i = 0; i!= u1bound.rows(); ++i) 
        {
            DenseVector<ResPrec_> actual_row = u1bound[i];
            for(typename DenseVector<ResPrec_>::ElementIterator j(actual_row.begin_elements()),
                                                            j_END(actual_row.end_elements());
                                                            //k2((*_v).begin_elements());
                                                                j!= j_END; ++j)
            {
                DenseVector<ResPrec_> flow =_flow_x<ResPrec_>(i,j.index());

                (*_v)[k2.index()] = flow[0];
                (*_v)[(k2.index())+1] = flow[1];
                (*_v)[(k2.index())+2] = flow[2];
                ++k2; ++k2; ++k2;
            }
        }
        */
       
        DenseVector<ResPrec_> *uFlow = _u->copy();
        _flow_x(*uFlow);
        _v = uFlow;
        cout << "v^T after building:\n";
        cout << stringify(*_v) << endl;
        /*OBSOLETE
        typename DenseVector<ResPrec_>::ElementIterator k3(_w->begin_elements());
        for (ulint i = 0; i!= u2bound.rows(); ++i) 
        {
            DenseVector<ResPrec_> actual_row = u2bound[i];
            for(typename DenseVector<ResPrec_>::ElementIterator j(actual_row.begin_elements()),
                                                            j_END(actual_row.end_elements());
                                                            //k3((*_w).begin_elements());
                                                                j!= j_END; ++j)
            {
                DenseVector<ResPrec_> flow = this->_flow_y<ResPrec_>(i,j.index());

                (*_w)[k3.index()] = flow[0];
                (*_w)[(k3.index())+1] = flow[1];
                (*_w)[(k3.index())+2] = flow[2];
                ++k3; ++k3; ++k3;
            }
        }
        */

        DenseVector<ResPrec_>* u2Flow = _u->copy();
        _flow_y(*u2Flow);
        _w = u2Flow;
        cout << "w^T after building:\n";
        cout << stringify(*_w) << endl;
       
  
        ///Now, that the relaxation vectors have been provided, the only thing left to do is to 
        ///compute the bottom slopes.
        typename DenseVector<ResPrec_>::ElementIterator l(_bottom_slopes_x->begin_elements());
        typename DenseVector<ResPrec_>::ElementIterator k4(_bottom_slopes_y->begin_elements());
        for (ulint i = 0; i!= bbound.rows(); ++i) 
        {
            DenseVector<ResPrec_> actual_row = bbound[i];
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
    
        cout << "Slopes after building:\n";
        cout << stringify(*_bottom_slopes_x) << endl;
        cout << stringify(*_bottom_slopes_y) << endl;
        std::cout << "Finished preprocessing.\n";
        
        //_do_postprocessing(1);

    }   

///Implementation of flow-processing functions.

    /**
     *
     * Flow computation in x-direction.
     * \param vector The vector, for which the flow is going to be computed.
     **/

    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>
    template <typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_flow_x(DenseVector<WorkPrec_> & vector)
    {
	if (!(vector.size() % 3)) 
	{
	    typename DenseVector<WorkPrec_>::ElementIterator resultvectoriterator(vector.begin_elements());
	    WorkPrec_ resultcomponentone, resultcomponenttwo, resultcomponentthree, gravterm;
	    for (typename DenseVector<WorkPrec_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()); l != l_end; ++l)
	    {	
	        // Compute additional gravitation-based term for flowcomponent two
	        gravterm = WorkPrec_(9.81 * (*l) * (*l) * 0.5);

	        // Compute the influence of the waterdepth
                if(*l!=0)
                {
	            resultcomponenttwo = 1 / *l;
	            resultcomponentthree = 1 / *l;
                }
                else
                {
	            resultcomponenttwo = 0;
	            resultcomponentthree = 0;
        
                }
	        ++l;

	        // Compute the influence of the waterflow in X-direction
		resultcomponentone = *l;
	        resultcomponenttwo = resultcomponenttwo * (*l) * (*l) + gravterm;
		resultcomponentthree *= *l;
	        ++l;

	        // Compute the influence of the waterflow in Y-direction and add the gravition-based term
	        resultcomponentthree *= *l ;

	        // Write the computed values into the resultvector
	        *resultvectoriterator = resultcomponentone;
	        ++resultvectoriterator;
	        *resultvectoriterator = resultcomponenttwo;
	        ++resultvectoriterator;
	        *resultvectoriterator = resultcomponentthree;
	        ++resultvectoriterator;
	    }
	}
	else
	{
	    std::cout << "Tststs... size of given vector does not match the requirements (size modulo 3 = 0).";
	}
        std::cout << "Finished Flow.\n";
    }


    /**
     *
     * Flow copmutation in y-direction.
     *
     * \param vector The vector, for which the flow is going to be computed.
     **/
  
    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>
    template <typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_flow_y(DenseVector<WorkPrec_> & vector)
    {
	if (!(vector.size() % 3)) 
	{
	    typename DenseVector<WorkPrec_>::ElementIterator resultvectoriterator(vector.begin_elements());
	    WorkPrec_ resultcomponentone, resultcomponenttwo, resultcomponentthree, gravterm;
	    for (typename DenseVector<WorkPrec_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()); l != l_end; ++l)
	    {
	        // Initialize locale resultcomponent variables
	        resultcomponentone = WorkPrec_(1);
	        resultcomponenttwo = WorkPrec_(1);
	        resultcomponentthree = WorkPrec_(1);
	
	        // Compute additional gravitation-based term for flowcomponent two
	        gravterm = WorkPrec_(9.81 * (*l) * (*l) * 0.5);

	        // Compute the influence of the waterdepth
                if(*l!=0)
                {
	            resultcomponenttwo = 1 / *l;
	            resultcomponentthree = 1 / *l;
                }
                else
                {
 	            resultcomponenttwo = 0;
	            resultcomponentthree = 0;
                           
                }
	        ++l;

	        // Compute the influence of the waterflow in X-direction
	        resultcomponenttwo *= *l;
	        ++l;

	        // Compute the influence of the waterflow in Y-direction and add the gravition-based term
	        resultcomponentone *= *l;
	        resultcomponenttwo *= *l;
	        resultcomponentthree = (resultcomponentthree * (*l) * (*l)) + gravterm;

	        // Write the computed values into the resultvector
	        *resultvectoriterator = resultcomponentone;
	        ++resultvectoriterator;
	        *resultvectoriterator = resultcomponenttwo;
	        ++resultvectoriterator;
	        *resultvectoriterator = resultcomponentthree;
	        ++resultvectoriterator;
                
	    }
	}
	else
	{
	    std::cout << "Tststs... size of given vector does not match the requirements (size modulo 3 = 0).";
	}
        std::cout << "Finished Flow.\n";
    }

    /**
     *
     * Flow computation in x-direction for a single cell.
     *
     * This function is used only by the preprocessing stage.
     *
     * \param i First coordinate of the processed cell.
     * \param j Second coordinate of the processed cell.
     *
     **/

    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>	  
    template <typename WorkPrec_>
    DenseVector<WorkPrec_> RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_flow_x(uint i, uint j)
    {
	DenseVector<WorkPrec_> result(ulint(3), ulint(0));
	WorkPrec_ temp = (*_v)[(_d_width + 4) * 3 * i + 3 * j];

        WorkPrec_ gravterm = 9.81 * temp * temp / 2;

	result[1] = 1 / temp;
        result[2] = result[1];

	temp = (*_v)[(_d_width + 4) * 3 * i + 3 * j + 1];
        result[0] = temp;
        result[1] *= temp * temp;
	result[2] *= temp;

        temp = (*_v)[(_d_width + 4) * 3 * i + 3 * j + 2];
	result[1] += gravterm;
        result[2] *= temp;
        std::cout << "Finished Flow x (cell). \n";

	return result;

    }

    /**
     *
     * Flow computation in y-direction for a single cell.
     *
     * This function is used only by the preprocessing stage
     *
     * \param i First coordinate of the processed cell.
     * \param j Second coordinate of the processed cell.
     *
     **/

    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>	  
    template <typename WorkPrec_>
    DenseVector<WorkPrec_> RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_flow_y(uint i, uint j)
    {
        DenseVector<WorkPrec_> result(ulint(3), ulint(0));
	WorkPrec_ temp = (*_w)[(_d_width + 4) * 3 * i + 3 * j];

	WorkPrec_ gravterm = 9.81 * temp * temp / 2;

        result[1] = 1 / temp;
	result[2] = result[1];

        temp = (*_w)[(_d_width + 4) * 3 * i + 3 * j + 1];
	result[1] *= temp;

        temp = (*_w)[(_d_width + 4) * 3 * i + 3 * j + 2];
	result[0] *= temp;
        result[1] *= temp;
        result[2] = result[2] * temp * temp + gravterm;
        std::cout << "Finished Flow y (cell).\n";

	return result;

    }


    /**
     *
     * Source term computation for a single cell (i, j).
     *
     * \param i First coordinate of the processed cell.
     * \param j Second coordinate of the processed cell.
     *
     **/

    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>
    template <typename WorkPrec_>
    DenseVector<WorkPrec_> RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_source(uint i, uint j)
    {
	DenseVector<WorkPrec_> result(ulint(3), ulint(0), ulint(1));
	WorkPrec_ h = (*_u)[(_d_width + 4) * 3 * i + 3 * j];
	if (h > 0)
	{
	    WorkPrec_ q1 = (*_u)[(_d_width + 4) * 3 * i + 3 * j + 1];
	    WorkPrec_ q2 = (*_u)[(_d_width + 4) * 3 * i + 3 * j + 2];
	    
	    result[0] = 0;
	    result[1] = _manning_n_squared * pow(h, -7/3) * sqrt(q1 * q1 + q2 * q2) * (-1);
	    result[2] = result[1];
    
	    result[1] = ((result[1] * q1) - (h * (*_bottom_slopes_x)[(_d_width + 4)* i + j])) * 9.81;
	    result[2] = ((result[2] * q2) - (h * (*_bottom_slopes_y)[(_d_width + 4)* i + j])) * 9.81;
            std::cout << "Finished simple flow.\n";
 
	    return result;
    	}
	else
	{
	    result[0] = 0;
	    result[1] = 0;
	    result[2] = 0;
	    return result;
	}

    }

    /** 
     *
     * Source term computation for a whole vector.
     *
     * \param vector Densevector for which the flow should be computed.
     *
     **/

    template <typename ResPrec_,
	      typename PredictionPrec1_,
	      typename PredictionPrec2_,
	      typename InitPrec1_,
	      typename InitPrec2_>
    template <typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_source(DenseVector<WorkPrec_>& vector)
    {
        WorkPrec_ oneThird = WorkPrec_(-1./3./*0.33333333333333333333333333333333333333333333333333333333333333333333333333333*/);
	if (!(vector.size() % 3)) 
	{
	    typename DenseVector<WorkPrec_>::ElementIterator resultvectoriterator(vector.begin_elements());
    	    typename DenseVector<WorkPrec_>::ConstElementIterator bottomslopesxiterator(_bottom_slopes_x->begin_elements());
	    typename DenseVector<WorkPrec_>::ConstElementIterator bottomslopesyiterator(_bottom_slopes_y->begin_elements());
	    WorkPrec_ h, u1, u2, friction;
	    for (typename DenseVector<WorkPrec_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()); l != l_end; ++l)
	    {
	        // Fetch values for source term computation
	        h = *l;
		++l;
                if(h!=0)
                {
	            u1 = *l/h ;
		    ++l;
	            u2 = *l/h ;
	        }
                else
                {
                    u1 = 0;
                    u2 = 0;
                    ++l;
                }
	        // Compute the influence of the waterdepth
	        *resultvectoriterator = WorkPrec_(0);
		++resultvectoriterator;
	        if(h>=0)
                {
		    friction = this->_manning_n_squared * pow(h, oneThird) * sqrt(u1 * u1 + u2 * u2) * (-1);
                    cout << "h:"<< stringify(h)<<"pow:" << stringify(  pow(h, oneThird)  ) << endl;
                }
                else
                {
                    friction = 0;
                }
                *resultvectoriterator = ((friction * u1) - (h * *bottomslopesxiterator)) * 9.81;
		++bottomslopesxiterator;
		++resultvectoriterator;
		
		*resultvectoriterator = ((friction * u2) - (h * *bottomslopesyiterator)) * 9.81;
		++bottomslopesyiterator;
		++resultvectoriterator;
	    }
	}
	else
	{
	    std::cout << "Tststs... size of given vector does not match the requirements (size modulo 3 = 0).";
	}
	std::cout << "Finished Source.\n";
    }

    /**
     *
     * Matrix assembly.
     *
     * \param  m1 The first matrix to assemble.
     * \param  m3 The second matrix to assemble.
     *
     **/
    template<typename ResPrec_ , typename PredictionPrec1_, typename PredictionPrec2_, typename InitPrec1_, typename InitPrec2_>
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>
        ::_assemble_matrix1(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& m3, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v)
    {
        ///The bands containing data.
        DenseVector<WorkPrec_> m1diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m1bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m1bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m1bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m3diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m3bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m3bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m3bandMinus1(_u->size(),ulint( 0));//three
        m1.band(ulint(0)) = m1diag;
        m1.band(ulint(3)) = m1bandPlus1;
        m1.band(ulint(6)) = m1bandPlus2;
        m1.band(ulint(-3)) = m1bandMinus1;
        m3.band(ulint(0)) = m3diag;
        m3.band(ulint(3)) = m3bandPlus1;
        m3.band(ulint(6)) = m3bandPlus2;
        m3.band(ulint(-3)) = m3bandMinus1;
 
        
        ///Necessary values to be temporarily saved.
        DenseVector<WorkPrec_> tempPlus(ulint(3),ulint(0));
        DenseVector<WorkPrec_> tempTopPlus(ulint(3),ulint(0));

        DenseVector<WorkPrec_> phiPlusOld(ulint(3),ulint(0));
        DenseVector<WorkPrec_> phiPlusNew(ulint(3),ulint(0));
        DenseVector<WorkPrec_> phiMinusNew(ulint(3),ulint(0));

        DenseVector<WorkPrec_> tempMinus(ulint(3),ulint(0));
        DenseVector<WorkPrec_> tempTopMinus(ulint(3),ulint(0));
        
        WorkPrec_ phiMinusOld;
        WorkPrec_ temp;
        WorkPrec_ tempTop;
        WorkPrec_ prefac = _delta_t/4*_delta_x;

        ///Needed Iterators.
        typename DenseVector<WorkPrec_>::ElementIterator d(m1diag.begin_elements());
        typename DenseVector<WorkPrec_>::ConstElementIterator i(u->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b1(m1bandPlus1.begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b2(m1bandPlus2.begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator bminus1(m1bandMinus1.begin_elements());
        
        ///Iterate through the vectors in order to avoid boundary access.
        for( ; i.index() < 6*(_d_width+4); ++i);
        for( ; d.index() < 6*(_d_width+4); ++d);                               
        for( ; b1.index() < 6*(_d_width+4); ++b1);                               
        for( ; b2.index() < 6*(_d_width+4); ++b2);
        for( ; bminus1.index() < 6*(_d_width+4); ++bminus1);                               
        
        //while(i.index() < 3*(_d_width+4) * (_d_height))
        //{
            for(unsigned long k=0; k<3; ++k)
            {
                tempPlus[k]= (*v)[i.index()] + (*_c)[k]*(*u)[i.index()];
                /*cout << "u[i]:" << stringify((*u)[i.index()]) << endl;
                cout << "v[i]:" << stringify((*v)[i.index()]) << endl;
                cout << "temPlus:"<< stringify(tempPlus[k])<<endl;*/
                ++i;++d;++b1;++b2;++bminus1;
            }

            for(unsigned long k=0; k<3; ++k)
            {
                temp= (*v)[i.index()] + (*_c)[k]*(*u)[i.index()];
                tempTopPlus[k] = temp - tempPlus[k];
                tempPlus[k] = temp;
                tempMinus[k] =  (*v)[i.index()] - (*_c)[k]*(*u)[i.index()];
                ++i;++d;++b1;++b2;++bminus1;
            }

            for(unsigned long k=0; k<3; ++k)
            {
                temp= (*v)[i.index()] - (*_c)[k]*(*u)[i.index()];
                tempTopMinus[k] = temp - tempMinus[k];
                tempMinus[k] = temp;
                temp= (*v)[i.index()] + (*_c)[k]*(*u)[i.index()];
                tempTop = temp - tempPlus[k];
                if(tempTop != 0)
                {
                    phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                }
                else
                {
                    phiPlusOld[k] = WorkPrec_(0);
                }
                tempPlus[k]=temp;
                tempTopPlus[k]=tempTop;
                ++i;
            }

            for(unsigned long k=0; k<3; ++k)
            {
                temp = (*v)[i.index()] - (*_c)[k]*((*u)[i.index()]); //temp = v(i) - c(k)*u(i);
                tempTop = temp - tempMinus[k]; //temp_top = temp - temp_minus(k);
                if(tempTop != 0)
                {
                    phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);//phi_minus_new(k) = Limiter(temp_top_minus(k) / temp_top);
                }
                else
                {
                    phiMinusNew[k] = WorkPrec_(0);
                }
                tempMinus[k]=temp;//switch(temp, temp_minus(k));
                tempTopMinus[k] = tempTop;//switch(temp_top, temp_top_minus(k));
                temp  = (*v)[i.index()] + (*_c)[k]* (*u)[i.index()];//temp = v(i) + c(k)*u(i);
                tempTop = temp - tempPlus[k];//temp_top = temp - temp_plus(k);
                if(tempTop != 0)
                {
                    phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);//phi_plus_new(k) = Limiter(temp_top_plus(k) / temp_top);
                }
                else
                {
                    phiPlusNew[k] = 0;
                }
                tempPlus[k]= temp;//switch(temp, temp_plus(k));
                tempTopPlus[k] = tempTop;//switch(temp_top, temp_top_plus(k));
                ++i;
            }

            for(unsigned long x = 0; x < _u->size(); ++x)
            {
                for(unsigned long k =0; k<3; ++k)
                {
                    temp = prefac *(2 - phiPlusOld[k]);
                   
                    m1bandMinus1[bminus1.index()] =temp;
                    m3bandMinus1[bminus1.index()] =temp * (*_c)[k];
        
                    m1diag[d.index()] = prefac * (phiPlusNew[k] + phiPlusOld[k] + phiMinusNew[k]);
                    m3diag[d.index()] = (*_c)[k] * prefac *(4 - phiPlusNew[k] - phiPlusOld[k] + phiMinusNew[k]);
                    
                    phiPlusOld[k]= phiPlusNew[k];
                    phiMinusOld = phiMinusNew[k];
                    temp = (*v)[i.index()] - (*_c)[k]*(*u)[i.index()]; //temp = v(i) - c(k)*u(i);
                    tempTop = temp - tempMinus[k]; //temp_top = temp - temp_minus(k);
                    
                    if(tempTop != 0)
                    {
                        phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);//phi_minus_new(k) = Limiter(temp_top_minus(k) / temp_top);
                    }
                    else
                    {
                        phiMinusNew[k] = WorkPrec_(0);
                    }
                    
                    tempMinus[k]=temp;//switch(temp, temp_minus(k));
                    tempTopMinus[k] = tempTop;//switch(temp_top, temp_top_minus(k));
                    temp  = (*v)[i.index()] + (*_c)[k]* (*u)[i.index()];//temp = v(i) + c(k)*u(i);
                    tempTop = temp - tempPlus[k];//temp_top = temp - temp_plus(k);
                    if(tempTop != 0)
                    {
                        phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);//phi_plus_new(k) = Limiter(temp_top_plus(k) / temp_top);
                    }
                    else
                    {
                        phiPlusNew[k] = 0;
                    }
                    
                    tempPlus[k]= temp;//switch(temp, temp_plus(k));
                    tempTopPlus[k] = tempTop;//switch(temp_top, temp_top_plus(k));
                    m1bandPlus1[b1.index()] = prefac * (-2 - phiPlusOld[k] - phiMinusOld - phiMinusNew[k]);
                    m3bandPlus1[b1.index()] = (*_c)[k]*prefac * (2 - phiPlusOld[k] + phiMinusOld + phiMinusNew[k]);
                    m1bandPlus2[b2.index()] = prefac* phiMinusNew[k];
                    m3bandPlus2[b2.index()] = (*_c)[k] * prefac *(-phiMinusNew[k]);
                    ++i;//++d;++b1;++b2;++bminus1;
                }
                
 
                //}
                ++d;++b1;++b2;++bminus1;
                ++d;++b1;++b2;++bminus1;
                ++d;++b1;++b2;++bminus1;
                ++d;++b1;++b2;++bminus1;
                ++d;++b1;++b2;++bminus1;
                ++d;++b1;++b2;++bminus1;
                
            }

 /*           m1.band(ulint(0)) = m1diag;
            m1.band(ulint(3)) = m1bandPlus1;
            m1.band(ulint(6)) = m1bandPlus2;
            m1.band(ulint(-3)) = m1bandMinus1;
            m3.band(ulint(0)) = m3diag;
            m3.band(ulint(3)) = m3bandPlus1;
            m3.band(ulint(6)) = m3bandPlus2;
            m3.band(ulint(-3)) = m3bandMinus1;*/
            
            cout << "M_1:" << stringify(m1.band(ulint(0))) << endl;
            cout << "M_1:" << stringify(m1.band(ulint(3))) << endl;
            cout << "M_1:" << stringify(m1.band(ulint(6))) << endl;
            cout << "M_1:" << stringify(m1.band(ulint(-3))) << endl;
            std::cout << "Finished Matrix Assembly 1.\n";
        }
    /**
      * First setup of values.
      * 
      *
      **/
    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>:: _do_setup_stage1(DenseVector<WorkPrec_>& su, DenseVector<WorkPrec_>& sv, DenseVector<WorkPrec_>& sw)
    {
       
        /*_u_temp = new DenseVector<WorkPrec_>(*(_u->copy()));//reinterpret_cast<DenseVector<WorkPrec_>*>(_u);
        _v_temp = new DenseVector<WorkPrec_>(*(_v->copy()));//reinterpret_cast<DenseVector<WorkPrec_>*>(_v);
        _w_temp = new DenseVector<WorkPrec_>(*(_w->copy()));//reinterpret_cast<DenseVector<WorkPrec_>*>(_w);
        */
        WorkPrec_ prefac;
        if(_eps != _delta_t)
        {
            prefac = 1/(_eps - _delta_t);
        }
        else
        {
	  cout << "prefac is invalid!\n";
        }
        DenseVector<WorkPrec_> *v = _v->copy();
        DenseVector<WorkPrec_> *u1 = _u->copy();
        DenseVector<WorkPrec_> *u2 = _u->copy();
        //DenseVector<WorkPrec_> flow1(_u->size(),ulint(0));
        _flow_x<WorkPrec_>(*u1);
        DenseVector<WorkPrec_> tempsum = VectorScaledSum<>::value(*v, *u1, _eps,_delta_t);
        delete u1;
        sv = ScalarVectorProduct<>::value(prefac,tempsum);
        DenseVector<WorkPrec_> *w = _w->copy();
        //DenseVector<WorkPrec_> flow2(_u->size(), ulint(0), ulint(1));
        _flow_y<WorkPrec_>(*u2);
        DenseVector<WorkPrec_> tempsum2 = VectorScaledSum<>::value(*w, *u2, _eps,_delta_t);
        delete u2;
        
        sw = ScalarVectorProduct<>::value(prefac,tempsum2);
        
        cout << "Temp relax vectors after building:\n";
        cout << stringify(*_u_temp) << endl;
        cout << stringify(*_v_temp) << endl;
        cout << stringify(*_w_temp) << endl;
        std::cout << "Finished Setup 1.\n";
    }

    /** The prediction stage.
      *
      * \param predictedu Temp u.
      * \param predictedv Temp v.
      * \param predictedw Temp w.
      **/
    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    template<typename WorkPrec_>
        void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>:: _do_prediction(DenseVector<WorkPrec_>& predictedu, DenseVector<WorkPrec_>& predictedv, DenseVector<WorkPrec_>& predictedw)
    {
        BandedMatrix<WorkPrec_> m1(_u->size());
        BandedMatrix<WorkPrec_> m2(_u->size());
        BandedMatrix<WorkPrec_> m3(_u->size());
        BandedMatrix<WorkPrec_> m4(_u->size());
        
        
        /*BandedMatrix<WorkPrec_>* m1= new BandedMatrix<WorkPrec_>(_u->size());
        BandedMatrix<WorkPrec_>* m2= new BandedMatrix<WorkPrec_>(_u->size());
        BandedMatrix<WorkPrec_>* m3= new BandedMatrix<WorkPrec_>(_u->size());
        BandedMatrix<WorkPrec_>* m4= new BandedMatrix<WorkPrec_>(_u->size());
        */
        _assemble_matrix1_DEBUG<WorkPrec_>(m1, m3, &predictedu, &predictedv);
        _assemble_matrix2_DEBUG<WorkPrec_>(m2, m4, &predictedu, &predictedw);
        
        BandedMatrix<WorkPrec_> m5 = *(m3.copy());
         
        BandedMatrix<WorkPrec_> m6(_u->size());
        _quick_assemble_matrix2<WorkPrec_>(m1, m6);
 
        BandedMatrix<WorkPrec_> m7 = *(m4.copy());
      
        BandedMatrix<WorkPrec_> m8(_u->size());
        _quick_assemble_matrix4<WorkPrec_>(m2, m8); 
	cout << "Prediction: Finished assembly.\n";
        //BandedMatrix<WorkPrec_> m5 = _quick_assemble_matrix1<WorkPrec_>(m3);
        //BandedMatrix<WorkPrec_> m6 = _quick_assemble_matrix2<WorkPrec_>(m1);
        //BandedMatrix<WorkPrec_> m7 = _quick_assemble_matrix3<WorkPrec_>(m4);
        //BandedMatrix<WorkPrec_> m8 = _quick_assemble_matrix4<WorkPrec_>(m2);
        
        DenseVector<WorkPrec_>* tempu = predictedu.copy();
        DenseVector<WorkPrec_>* tempv = predictedv.copy();
        DenseVector<WorkPrec_>* tempw = predictedw.copy();
	cout << "Prediction: Before matrix*vector.\n";
        DenseVector<WorkPrec_> temp1 = MatrixVectorProduct<>::value<WorkPrec_,WorkPrec_>(m1,*tempv);
        DenseVector<WorkPrec_> *tempu2 = predictedu.copy();
	cout << "First product solved.\n";
        DenseVector<WorkPrec_> temp2 = MatrixVectorProduct<>::value<WorkPrec_,WorkPrec_>(m2,*tempw);
	cout << "Second product solved.\n";
	
        DenseVector<WorkPrec_> temp3 = MatrixVectorProduct<>::value<WorkPrec_,WorkPrec_>(m3,*tempu);
	cout << "Third product solved.\n";
        DenseVector<WorkPrec_> temp4 = MatrixVectorProduct<>::value<WorkPrec_,WorkPrec_>(m4,*tempu2);
        cout << "Fourth product solved.\n";
	DenseVector<WorkPrec_>* source = predictedu.copy(); 
        _source(*source);
	cout << "Source solved.\n";
	DenseVector<WorkPrec_>* predicteduTemp = predictedu.copy();
	*predicteduTemp = VectorSum<>::value<WorkPrec_,WorkPrec_>(*predicteduTemp, temp1);
        *predicteduTemp = VectorSum<>::value<WorkPrec_,WorkPrec_>(*predicteduTemp,temp2);
        *predicteduTemp = VectorSum<>::value(*predicteduTemp, temp3);
        *predicteduTemp = VectorSum<>::value(*predicteduTemp, temp4);
        *predicteduTemp = VectorSum<>::value(*predicteduTemp, *source);
        delete source;
	
	cout << "First accu solved.\n";
        DenseVector<WorkPrec_> *tempu3 = predictedu.copy();
        DenseVector<WorkPrec_> *tempv2 = predictedv.copy();
        DenseVector<WorkPrec_> *tempw2 = predictedw.copy();
        
        
        DenseVector<WorkPrec_> temp11 = MatrixVectorProduct<>::value(m5,*tempv2);
        DenseVector<WorkPrec_> *tempu4 = predictedu.copy();
        DenseVector<WorkPrec_> temp22 = MatrixVectorProduct<>::value(m6,*tempu3);
        DenseVector<WorkPrec_> temp33 = MatrixVectorProduct<>::value(m7,*tempw2);
        DenseVector<WorkPrec_> temp44 = MatrixVectorProduct<>::value(m8,*tempu4);
        
	predictedu = *predicteduTemp;
	predictedv = VectorSum<>::value(*(predictedv.copy()), temp11);
        predictedv = VectorSum<>::value(*(predictedv.copy()), temp22);
        predictedw = VectorSum<>::value(*(predictedw.copy()), temp33);
        predictedw = VectorSum<>::value(*(predictedw.copy()), temp44);
        
        /*delete m1;
        delete m2;
        delete m3;
        delete m4;
        delete m5;
        delete m6;
        delete m7;
        delete m8;*/
        std::cout << "Finished Prediction.\n";
        
    }

    /** Second setup of values.
      *
      * \param predictedu Temp u.
      * \param predictedv Temp v.
      * \param predictedw Temp w.
      **/
    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>:: _do_setup_stage2(DenseVector<WorkPrec_>& predictedu, DenseVector<WorkPrec_>& predictedv, DenseVector<WorkPrec_>& predictedw)
    {
        ///Apply flow to newest u:
        DenseVector<WorkPrec_>* f = predictedu.copy();
        _flow_x(*f);

        ///Compute linear combinations und accumulate:
        DenseVector<WorkPrec_> innersum1 = VectorScaledSum<>::value<WorkPrec_, WorkPrec_, double, WorkPrec_>(predictedv,*f, _eps, _delta_t);
        delete f;
        ///Apply flow to old u:
        DenseVector<WorkPrec_>* flow = _u_temp->copy();
        _flow_x(*flow);

        DenseVector<WorkPrec_> innersum2 = VectorScaledSum<>::value<WorkPrec_, WorkPrec_, WorkPrec_, WorkPrec_>(*_v_temp,*flow, -2*_delta_t, 2*_delta_t);
        delete flow;
        predictedv = VectorSum<>::value<WorkPrec_, WorkPrec_>(innersum1, innersum2);
        ///Scale sum:
        predictedv = ScalarVectorProduct<>::value(1+(_eps/_delta_t), predictedv);
        
        ///Repeat for w:
        DenseVector<WorkPrec_>* flow2 = predictedu.copy();
        _flow_y(*flow2);
    
        innersum1 = VectorScaledSum<>::value<WorkPrec_, WorkPrec_, double, WorkPrec_>(predictedw, *flow2, _eps, _delta_t);
        delete flow2;
        DenseVector<WorkPrec_>* flow3 = _u_temp->copy();
        _flow_x(*flow3);

        innersum2 = VectorScaledSum<>::value<WorkPrec_, WorkPrec_, WorkPrec_, WorkPrec_>(*_w_temp, *flow3, -2*_delta_t, 2*_delta_t);
        delete flow3;
        predictedw = VectorSum<>::value<WorkPrec_, WorkPrec_>(innersum1, innersum2);
        predictedw = ScalarVectorProduct<>::value(1+(1/_delta_t), predictedw);
        std::cout << "Finished Setup 2.\n";
    } 
    
    /** Implementation of the correction stage.
      *
      *\param predictedu Temp u.
      *\param predictedv Temp v.
      *\param predictedw Temp w.
      *
      **/
    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>:: _do_correction(DenseVector<ResPrec_>& predictedu,
    DenseVector<ResPrec_>& predictedv,
    DenseVector<ResPrec_>& predictedw)
    {   
        ///correct first 2(w+4)+2 ghost cells (tripels)
        typename DenseVector<ResPrec_>::ElementIterator iter(_u->begin_elements());
        while(iter.index()<(6*(_d_width+4)+6))
        {
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
            (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
            (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);


        }
        cout << stringify(iter.index()) << endl;
        
        unsigned long count =0;//if made w steps, ignore four.
        ///Iterate through predicted u,v,w - vectors, compute weighted sum , read out h_ij, care about ghost cells.
        typename DenseMatrix<ResPrec_>::ElementIterator h(_height->begin_elements());
        while(iter.index()<3*((_d_height+2)*(_d_width+4)))
        {
            /*PredictionPrec2_ a = (*_u)[iter.index()];
            (*_u)[iter.index()] = 0.5*(predictedu[iter.index()] + a);
            cout << stringify(predictedu[iter.index()])<< "+" << stringify(a) << "/2 =" << stringify((*_u)[iter.index()] )<<endl;*/
            if(count < _d_width)
            {
                PredictionPrec2_ a = (*_u)[iter.index()];
                (*_u)[iter.index()] = 0.5*(predictedu[iter.index()] + a);
                cout << stringify(predictedu[iter.index()])<< "+" << stringify(a) << "/2 =" << stringify((*_u)[iter.index()] )<<endl;
                
                *h = (*_u)[iter.index()];
                ++h;
                ++count;
              
                (*_v)[iter.index()] = 0.5*(predictedv[iter.index()]+ (*_v)[iter.index()]);
                (*_w)[iter.index()] = 0.5*(predictedw[iter.index()]+ (*_w)[iter.index()]);
                ++iter;
                (*_u)[iter.index()] = 0.5*(predictedu[iter.index()]+ (*_u)[iter.index()]);
                (*_v)[iter.index()] = 0.5*(predictedv[iter.index()]+ (*_v)[iter.index()]);
                (*_w)[iter.index()] = 0.5*(predictedw[iter.index()]+ (*_w)[iter.index()]);
                ++iter;
                (*_u)[iter.index()] = 0.5*(predictedu[iter.index()]+ (*_u)[iter.index()]);
                (*_v)[iter.index()] = 0.5*(predictedv[iter.index()]+ (*_v)[iter.index()]);
                (*_w)[iter.index()] = 0.5*(predictedw[iter.index()]+ (*_w)[iter.index()]);
                ++iter;
 
            }
            else
            {
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
                (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
                (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
                (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_u)[iter.index()] = 5;
                ++iter;
                (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
                (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
                (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
                count = 0;
            }
            cout << stringify(count)<<endl;
            
            
        }
        typename DenseVector<ResPrec_>::ElementIterator iter_END(_u->end_elements());
        while(iter!=iter_END)
        {
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_u)[iter.index()] = 5;
            ++iter;
            (*_v)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2];
            (*_v)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-2]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);
            (*_v)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-3] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-2] = (*_u)[iter.index()-3]*(*_u)[iter.index()-2]*(*_u)[iter.index()-1];
            (*_w)[iter.index()-1] = (*_u)[iter.index()-3]*(*_u)[iter.index()-1]*(*_u)[iter.index()-1]+ 
                                                (0.5*9.81*(*_u)[iter.index()-3]*(*_u)[iter.index()-3]);


           
        }
        std::cout << "Finished Correction.\n";
    }

    /** Capsule for the solution- computation in one timestep.
      *
      *
      **/
    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>:: solve()
    {
        _u_temp = _u->copy();
        _v_temp = _v->copy();
        _w_temp = _w->copy();
        
        /*_u_predicted = _u_temp->copy();
        _v_predicted = _v_temp->copy();
        _w_predicted = _w_temp->copy();
        */

        _do_setup_stage1<InitPrec1_>(*_u_temp, *_v_temp, *_w_temp );
        //DenseVector<PredictionPrec1_>* tu, * tv, * tw;
        /*tu = _u->copy();
        tv = _v->copy();
        tw = _w->copy();*/
        /*tu = new DenseVector<PredictionPrec1_>(*(_u_temp->copy()));
        tv = new DenseVector<PredictionPrec1_>(*(_v_temp->copy()));
        tw = new DenseVector<PredictionPrec1_>(*(_w_temp->copy()));
        */
        DenseVector<PredictionPrec1_>* predictedu = _u_temp->copy();
        DenseVector<PredictionPrec1_>* predictedv = _v_temp->copy();
        DenseVector<PredictionPrec1_>* predictedw = _w_temp->copy();
        
        _do_prediction<PredictionPrec1_>(*predictedu, *predictedv, *predictedw);
        //DenseVector<InitPrec2_> predictedu2(_u->size(), 0, 1);
        //DenseVector<InitPrec2_> predictedv2(_u->size(), 0, 1);
        //DenseVector<InitPrec2_> predictedw2(_u->size(), 0, 1);
        _do_setup_stage2<InitPrec2_>(*predictedu, *predictedv, *predictedw);
        _do_prediction<PredictionPrec2_>(*predictedu, *predictedv, *predictedw);
        cout << "Predicted u:\n";
        cout << stringify(*predictedu)<< endl;
        cout << "u before correction:\n";
        cout << stringify(*_u)<<endl;
        
        _do_correction(*predictedu, *predictedv, *predictedw);
        ++_solve_time;

        delete predictedu;
        delete predictedv;
        delete predictedw;
        delete _u_temp;
        delete _v_temp;
        delete _w_temp;
        /*predictedu = 0;
        predictedv = 0;
        predictedw = 0;
        _u_temp = 0;
        _v_temp = 0;
        _w_temp = 0;
        */
        cout << "Corrected u, finished solution, timestep:" << stringify(_solve_time) << endl;
        cout << stringify(*_u)<<endl;
    }

    /**
     *
     * Matrix assembly.
     *
     * \param  m2 The first matrix to assemble.
     * \param  m4 The second matrix to assemble.
     *
     **/
    template<typename ResPrec_ , typename PredictionPrec1_, typename PredictionPrec2_, typename InitPrec1_, typename InitPrec2_>
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>
        ::_assemble_matrix2(BandedMatrix<WorkPrec_>& m2, BandedMatrix<WorkPrec_>& m4, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* w)
    {   
        
        ///The bands containing data.
        DenseVector<WorkPrec_> m2diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m2bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m2bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m2bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m4diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m4bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m4bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m4bandMinus1(_u->size(),ulint( 0));//three
        
        ///Necessary values to be temporarily saved.
        DenseVector<WorkPrec_> tempPlus(ulint(3),ulint(0));
        DenseVector<WorkPrec_> tempTopPlus(ulint(3),ulint(0));

        DenseVector<WorkPrec_> phiPlusOld(ulint(3),ulint(0));
        DenseVector<WorkPrec_> phiPlusNew(ulint(3),ulint(0));
        DenseVector<WorkPrec_> phiMinusNew(ulint(3),ulint(0));

        DenseVector<WorkPrec_> tempMinus(ulint(3),ulint(0));
        DenseVector<WorkPrec_> tempTopMinus(ulint(3), ulint(0));
        
        WorkPrec_ phiMinusOld;
        WorkPrec_ temp;
        WorkPrec_ tempTop;
        WorkPrec_ prefac = _delta_t/4*_delta_y;

	for(unsigned long s=0; s< _d_width; ++s)
	{

	    ///Needed Iterators.
            typename DenseVector<WorkPrec_>::ElementIterator d(m2diag.begin_elements());
	    typename DenseVector<WorkPrec_>::ConstElementIterator i(u->begin_elements());
            typename DenseVector<WorkPrec_>::ElementIterator b1(m2bandPlus1.begin_elements());
	    typename DenseVector<WorkPrec_>::ElementIterator b2(m2bandPlus2.begin_elements());
	    typename DenseVector<WorkPrec_>::ElementIterator bminus1(m2bandMinus1.begin_elements());
            
            ///Iterate to the next column
            for(unsigned long f=0; f < 3*(s+2); ++f)
	    {
		++i;++d;++b1;++b2;++bminus1;
	    }

            
	    
	    for(unsigned long k=0; k<3; ++k)
            {
                tempPlus[k]= (*w)[i.index()] + (*_d)[k]*(*u)[i.index()];
		++i;++d;++b1;++b2;++bminus1;
	    }	

	    //Iterate to next column-element
	    for(unsigned long f=0; f<3*(_d_width+3);++f)
	    {
		++i;++d;++b1;++b2;++bminus1;
	    }

            for(unsigned long k=0; k<3; ++k)
            {
                temp= (*w)[i.index()] + (*_d)[k]*(*u)[i.index()];
                tempTopPlus[k] = temp - tempPlus[k];
                tempPlus[k] = temp;
                tempMinus[k] =  (*w)[i.index()] - (*_d)[k]*(*u)[i.index()];
                ++i;++d;++b1;++b2;++bminus1;
            }

	    //Iterate i to next column-element
	    for(unsigned long f=0; f<3*(_d_width+3);++f)
	    {
		++i;
	    }

            for(unsigned long k=0; k<3; ++k)
            {
                temp= (*w)[i.index()] - (*_d)[k]*(*u)[i.index()];
                tempTopMinus[k] = temp - tempMinus[k];
                tempMinus[k] = temp;
                temp= (*w)[i.index()] + (*_d)[k]*(*u)[i.index()];
                tempTop = temp - tempPlus[k];
                if(tempTop != 0)
                {
                    phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                }
                else
                {
                    phiPlusOld[k] = WorkPrec_(0);
                }
                tempPlus[k]=temp;
                tempTopPlus[k]=tempTop;
                ++i;
            }

	    //Iterate i to next column-element
	    for(unsigned long f=0; f<3*(_d_width+3);++f)
	    {
		++i;
	    }

            for(unsigned long k=0; k<3; ++k)
            {
                temp = (*w)[i.index()] - (*_d)[k]*((*u)[i.index()]); //temp = v(i) - c(k)*u(i);
                tempTop = temp - tempMinus[k]; //temp_top = temp - temp_minus(k);
                if(tempTop != 0)
                {
                    phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);//phi_minus_new(k) = Limiter(temp_top_minus(k) / temp_top);
                }
                else
                {
                    phiMinusNew[k] = WorkPrec_(0);
                }
                tempMinus[k]=temp;//switch(temp, temp_minus(k));
                tempTopMinus[k] = tempTop;//switch(temp_top, temp_top_minus(k));
                temp  = (*w)[i.index()] + (*_d)[k]* (*u)[i.index()];//temp = v(i) + c(k)*u(i);
                tempTop = temp - tempPlus[k];//temp_top = temp - temp_plus(k);
                if(tempTop != 0)
                {
                    phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);//phi_plus_new(k) = Limiter(temp_top_plus(k) / temp_top);
                }
                else
                {
                    phiPlusNew[k] = 0;
                }
                tempPlus[k]= temp;//switch(temp, temp_plus(k));
                tempTopPlus[k] = tempTop;//switch(temp_top, temp_top_plus(k));
                ++i;//++i;
            }
    
    
            for(unsigned long x = 0; x < _d_height; ++x)
            {

		//Iterate to next column-elements
		for(unsigned long f=0; f<3*(_d_width+3);++f)
		{
		    ++i;++d;++b1;++b2;++bminus1;
		}

                for(unsigned long k =0; k<3; ++k)
                {
                    temp = prefac *(2 - phiPlusOld[k]);
                       
                    m2bandMinus1[bminus1.index()] =temp;
                    m4bandMinus1[bminus1.index()] =temp * (*_d)[k];
          
                    m2diag[d.index()] = prefac * (phiPlusNew[k] + phiPlusOld[k] + phiMinusNew[k]);
                    m4diag[d.index()] = (*_d)[k] * prefac *(4 - phiPlusNew[k] - phiPlusOld[k] + phiMinusNew[k]);
                        
                    phiPlusOld[k]= phiPlusNew[k];
                    phiMinusOld = phiMinusNew[k];
                    temp = (*w)[i.index()] - (*_d)[k]*(*u)[i.index()]; //temp = v(i) - c(k)*u(i);
                    tempTop = temp - tempMinus[k]; //temp_top = temp - temp_minus(k);
                    if(tempTop != 0)
                    {
                        phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);//phi_minus_new(k) = Limiter(temp_top_minus(k) / temp_top);
                    }
                    else
                    {
                        phiMinusNew[k] = WorkPrec_(0);
                    }
                    tempMinus[k]=temp;//switch(temp, temp_minus(k));
                    tempTopMinus[k] = tempTop;//switch(temp_top, temp_top_minus(k));
                    temp  = (*w)[i.index()] + (*_d)[k]* (*u)[i.index()];//temp = v(i) + c(k)*u(i);
                    tempTop = temp - tempPlus[k];//temp_top = temp - temp_plus(k);
                    if(tempTop != 0)
                    {
                        phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);//phi_plus_new(k) = Limiter(temp_top_plus(k) / temp_top);
                    }
                    else
                    {
                        phiPlusNew[k] = 0;
                    }
                    tempPlus[k]= temp;//switch(temp, temp_plus(k));
                    tempTopPlus[k] = tempTop;//switch(temp_top, temp_top_plus(k));
                    m2bandPlus1[b1.index()] = prefac * (-2 - phiPlusOld[k] - phiMinusOld - phiMinusNew[k]);
                    m4bandPlus1[b1.index()] = (*_d)[k]*prefac * (2 - phiPlusOld[k] + phiMinusOld + phiMinusNew[k]);
                    m2bandPlus2[b2.index()] = prefac* phiMinusNew[k];
                    m4bandPlus2[b2.index()] = (*_d)[k] * prefac *(-phiMinusNew[k]);
                    ++i;++d;++b1;++b2;++bminus1;
                }  
            }
	}
        m2.band(ulint(0)) = m2diag;
        m2.band(ulint(3*(_d_width+4))) = m2bandPlus1;
        m2.band(ulint(6*(_d_width+4))) = m2bandPlus2;
        m2.band(ulint((-3)*(_d_width+4))) = m2bandMinus1;
        m4.band(ulint(0)) = m4diag;
        m4.band(ulint(3*(_d_width+4))) = m4bandPlus1;
        m4.band(ulint(6*(_d_width+4))) = m4bandPlus2;
        m4.band(ulint((-3)*(_d_width+4))) = m4bandMinus1;
        std::cout << "Finished Matrix Assembly 2.\n";
    }


    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    template<typename WorkPrec_>
    /*BandedMatrix<WorkPrec_>*/void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_quick_assemble_matrix2(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& result)
    {
	///Bands of the matrix which will be assembled.
	/*DenseVector<WorkPrec_> m6diag(_u->size(), ulint(0) ,ulint( 1));      //zero
        DenseVector<WorkPrec_> m6bandplus3(_u->size(), ulint(0) , ulint(1)); //one
        DenseVector<WorkPrec_> m6bandplus6(_u->size(), ulint(0) ,ulint( 1)); //two
        DenseVector<WorkPrec_> m6bandminus3(_u->size(), ulint(0) ,ulint(1));//three
        */
        /*DenseVector<WorkPrec_>* m6diag = new DenseVector<WorkPrec_>(_u->size(), ulint(0) ,ulint( 1));      //zero
        DenseVector<WorkPrec_>* m6bandplus3 = new DenseVector<WorkPrec_>(_u->size(), ulint(0) , ulint(1)); //one
        DenseVector<WorkPrec_>* m6bandplus6 = new DenseVector<WorkPrec_>(_u->size(), ulint(0) ,ulint( 1)); //two
        DenseVector<WorkPrec_>* m6bandminus3 = new DenseVector<WorkPrec_>(_u->size(), ulint(0) ,ulint(1));//three
        */

	DenseVector<WorkPrec_>* m6diag = (m1.band(ulint(0))).copy();
	DenseVector<WorkPrec_>* m6bandplus3 = (m1.band(ulint(3))).copy();
	DenseVector<WorkPrec_>* m6bandplus6 = (m1.band(ulint(6))).copy();
	DenseVector<WorkPrec_>* m6bandminus3 = (m1.band(ulint(-3))).copy();
        ///Needed Iterators.
        typename DenseVector<WorkPrec_>::ElementIterator d(m6diag->begin_elements());//, d_END(m6diag->end_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b1(m6bandplus3->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b2(m6bandplus6->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator bminus1(m6bandminus3->begin_elements());
        
	DenseVector<WorkPrec_> c_squared(*(_c->copy()));
	VectorElementwiseProduct<WorkPrec_>::value(c_squared, (*_c));

        for( ; d.index() < 6*(_d_width+4); ++d);                               
        for( ; b1.index() < 6*(_d_width+4); ++b1);                               
        for( ; b2.index() < 6*(_d_width+4); ++b2);
        for( ; bminus1.index() < 6*(_d_width+4); ++bminus1);                               

	while(d.index() < 3*(_d_width+4)*(_d_height-2))
	{
	    ++d; ++d; ++d; ++d; ++d; ++d;
	    ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
	    ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
	    ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;

	    for(ulint i = 0; i < _d_width; ++i)
	    {
		*d *= c_squared[0];
		*b1 *= c_squared[0];
		*b2 *= c_squared[0];
		*bminus1 *= c_squared[0];
		++d; ++b1; ++b2; ++bminus1;
		*d *= c_squared[1];
		*b1 *= c_squared[1];
		*b2 *= c_squared[1];
		*bminus1 *= c_squared[1];
		++d; ++b1; ++b2; bminus1;
		*d *= c_squared[2];
		*b1 *= c_squared[2];
		*b2 *= c_squared[2];
		*bminus1 *= c_squared[2];
		++d; ++b1; ++b2; ++bminus1;
	    }

	    ++d; ++d; ++d; ++d; ++d; ++d;
	    ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
	    ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
	    ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;
	}
	//BandedMatrix<WorkPrec_> result(m6diag.size());
	result.insert_band(0, m6diag);
	result.insert_band(3, m6bandplus3);
	result.insert_band(6, m6bandplus6);
	result.insert_band(-3, m6bandminus3);

	std::cout << "Finished Quick Assembly m2.\n";
 
	//return result;

    }

    template<typename ResPrec_,
             typename PredictionPrec1_,
             typename PredictionPrec2_,
             typename InitPrec1_,
             typename InitPrec2_> 
    template<typename WorkPrec_>
    /*BandedMatrix<WorkPrec_>*/void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_quick_assemble_matrix4(BandedMatrix<WorkPrec_>& m2, BandedMatrix<WorkPrec_>& result)
    {
	///Bands of the matrix which will be assembled.
	/*DenseVector<WorkPrec_> m8diag(_u->size(), ulint(0) ,ulint( 1));      //zero
        DenseVector<WorkPrec_> m8bandplus3(_u->size(), ulint(0) , ulint(1)); //one
        DenseVector<WorkPrec_> m8bandplus6(_u->size(), ulint(0) ,ulint( 1)); //two
        DenseVector<WorkPrec_> m8bandminus3(_u->size(), ulint(0) ,ulint(1));//three
        */

	DenseVector<WorkPrec_>* m8diag = (m2.band(ulint(0))).copy();
	DenseVector<WorkPrec_>* m8bandplus3 = (m2.band(ulint(3*(_d_width +4)))).copy();
	DenseVector<WorkPrec_>* m8bandplus6 = (m2.band(ulint(6*(_d_width +4)))).copy();
	DenseVector<WorkPrec_>* m8bandminus3 = (m2.band(ulint((-3)*(_d_width +4)))).copy();
        ///Needed Iterators.
        typename DenseVector<WorkPrec_>::ElementIterator d(m8diag->begin_elements()); //, d_END(m8diag->end_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b1(m8bandplus3->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator b2(m8bandplus6->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator bminus1(m8bandminus3->begin_elements());
 
	DenseVector<WorkPrec_> d_squared((*(_d->copy())));
	VectorElementwiseProduct<WorkPrec_>::value(d_squared, (*_d));

        for( ; d.index() < 6*(_d_width + 4); ++d);                               
        for( ; b1.index() < 6*(_d_width + 4); ++b1);                               
        for( ; b2.index() < 6*(_d_width + 4); ++b2);
        for( ; bminus1.index() < 6*(_d_width + 4); ++bminus1);                               

	while(d.index() < 3*(_d_width+4)*(_d_height-2))
	{
	    ++d; ++d; ++d; ++d; ++d; ++d;
	    ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
	    ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
	    ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;

	    for(ulint i = 0; i < _d_width; ++i)
	    {
		*d *= d_squared[0];
		*b1 *= d_squared[0];
		*b2 *= d_squared[0];
		*bminus1 *= d_squared[0];
		++d; ++b1; ++b2; ++bminus1;
		*d *= d_squared[1];
		*b1 *= d_squared[1];
		*b2 *= d_squared[1];
		*bminus1 *= d_squared[1];
		++d; ++b1; ++b2; bminus1;
		*d *= d_squared[2];
		*b1 *= d_squared[2];
		*b2 *= d_squared[2];
		*bminus1 *= d_squared[2];
		++d; ++b1; ++b2; ++bminus1;
	    }

	    ++d; ++d; ++d; ++d; ++d; ++d;
	    ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
	    ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
	    ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;
	}
	
	//BandedMatrix<WorkPrec_> result(m8diag.size());

	result.insert_band(0, m8diag);
	result.insert_band(3*(_d_width +4), m8bandplus3);
	result.insert_band(6*(_d_width +4), m8bandplus6);
	result.insert_band((-3)*(_d_width +4),m8bandminus3);
        std::cout << "Finished Quick Assembly m4.\n";	
 
	//return result;
	
    }
    /**
     *
     * Matrix assembly. DEBUG version.
     *
     * \param  m1 The first matrix to assemble.
     * \param  m3 The second matrix to assemble.
     *
     **/
    template<typename ResPrec_ , typename PredictionPrec1_, typename PredictionPrec2_, typename InitPrec1_, typename InitPrec2_>
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>
        ::_assemble_matrix1_DEBUG(BandedMatrix<WorkPrec_>& m1, BandedMatrix<WorkPrec_>& m3, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v)
    {
        ///The bands containing data.
        /*DenseVector<WorkPrec_> m1diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m1bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m1bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m1bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m3diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m3bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m3bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m3bandMinus1(_u->size(),ulint( 0));//three
        m1.band(ulint(0)) = m1diag;
        m1.band(ulint(3)) = m1bandPlus1;
        m1.band(ulint(6)) = m1bandPlus2;
        m1.band(ulint(-3)) = m1bandMinus1;
        m3.band(ulint(0)) = m3diag;
        m3.band(ulint(3)) = m3bandPlus1;
        m3.band(ulint(6)) = m3bandPlus2;
        m3.band(ulint(-3)) = m3bandMinus1;
        */
        DenseVector<WorkPrec_> m1diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m1bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m1bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m1bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m3diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m3bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m3bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m3bandMinus1(_u->size(),ulint( 0));//three

        typename DenseVector<WorkPrec_>::ElementIterator ui(u->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator ui_END(u->end_elements());
        typename DenseVector<WorkPrec_>::ElementIterator vi(u->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator vi_END(u->end_elements());
        
        //Preparing left operand for all scalar products
        DenseVector<WorkPrec_> constantVector(ulint(12), ulint(0));
        constantVector[0] = WorkPrec_(-1);
        constantVector[1] = WorkPrec_(-1);
        constantVector[2] = WorkPrec_(-1);
        constantVector[3] = WorkPrec_(1);
        constantVector[4] = WorkPrec_(1);
        constantVector[5] = WorkPrec_(1);
        constantVector[6] = WorkPrec_(-((*_c)[0]));
        constantVector[7] = WorkPrec_(-((*_c)[1]));
        constantVector[8] = WorkPrec_(-((*_c)[2]));
        constantVector[9]  = WorkPrec_((*_c)[0]);
        constantVector[10] = WorkPrec_((*_c)[1]);
        constantVector[11] = WorkPrec_((*_c)[2]);
        //Preparing left operand for all scalar products
        DenseVector<WorkPrec_> constantVectorMinus(ulint(12), ulint(0));
        constantVectorMinus[0] = WorkPrec_(-1);
        constantVectorMinus[1] = WorkPrec_(-1);
        constantVectorMinus[2] = WorkPrec_(-1);
        constantVectorMinus[3] = WorkPrec_(1);
        constantVectorMinus[4] = WorkPrec_(1);
        constantVectorMinus[5] = WorkPrec_(1);
        constantVectorMinus[6] = WorkPrec_((*_c)[0]);
        constantVectorMinus[7] = WorkPrec_((*_c)[1]);
        constantVectorMinus[8] = WorkPrec_((*_c)[2]);
        constantVectorMinus[9]  = WorkPrec_(-(*_c)[0]);
        constantVectorMinus[10] = WorkPrec_(-(*_c)[1]);
        constantVectorMinus[11] = WorkPrec_(-(*_c)[2]);

 
        while(ui!=ui_END)
        {
            //Prepare right operands for band_-1:
            //thetaXPlus_iMinus1_j
            DenseVector<WorkPrec_> rightMinus1Upper(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightMinus1Upper[0] = (*v)[vi.index()-6];
                rightMinus1Upper[1] = (*v)[vi.index()-5];
                rightMinus1Upper[2] = (*v)[vi.index()-4];
                rightMinus1Upper[3] = (*v)[vi.index()-3];
                rightMinus1Upper[4] = (*v)[vi.index()-2];
                rightMinus1Upper[5] = (*v)[vi.index()-1];

                rightMinus1Upper[6] = (*u)[ui.index()-6];
                rightMinus1Upper[7] = (*u)[ui.index()-5];
                rightMinus1Upper[8] = (*u)[ui.index()-4];
                rightMinus1Upper[9] = (*u)[ui.index()-3];
                rightMinus1Upper[10] = (*u)[ui.index()-2];
                rightMinus1Upper[11] = (*u)[ui.index()-1];
                



            }
            DenseVector<WorkPrec_> rightMinus1Lower(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6) 
            {
                rightMinus1Lower[0] = (*v)[vi.index()-3];
                rightMinus1Lower[1] = (*v)[vi.index()-2];
                rightMinus1Lower[2] = (*v)[vi.index()-1];
                rightMinus1Lower[3] = (*v)[vi.index()];
                rightMinus1Lower[4] = (*v)[vi.index()+1];
                rightMinus1Lower[5] = (*v)[vi.index()+2];

                rightMinus1Lower[6] = (*u)[ui.index()-3];
                rightMinus1Lower[7] = (*u)[ui.index()-2];
                rightMinus1Lower[8] = (*u)[ui.index()-1];
                rightMinus1Lower[9] = (*u)[ui.index()];
                rightMinus1Lower[10] = (*u)[ui.index()+1];
                rightMinus1Lower[11] = (*u)[ui.index()+2];
                
            }

            //Compute ScalarProducts and Theta - value for band_-1:
            WorkPrec_ minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightMinus1Upper, constantVector );
            WorkPrec_ minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightMinus1Lower, constantVector );
            WorkPrec_ thetaXPlus_iMinus1_j;
            if(minusLower!=0)
            {
                 thetaXPlus_iMinus1_j = minusUpper / minusLower;
            }
            //Compute limitation of Theta for band_-1:
            WorkPrec_ thetaXPlus_iMinus1_j_limited = min_mod_limiter(thetaXPlus_iMinus1_j);

            //Compute matrix element value:
            m1bandMinus1[ui.index()] = (_delta_t / 4*_delta_x) * (2 -  thetaXPlus_iMinus1_j_limited);
            m3bandMinus1[ui.index()] = (_delta_t / 4*_delta_x) * (2 -  thetaXPlus_iMinus1_j_limited);
            //FINISHED band_-1.

            //Prepare right operands for diagonal:
            //First Theta - value (thetaXPlus_i_j)
            DenseVector<WorkPrec_> rightDiagUpper(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagUpper[0] = (*v)[vi.index()-3];
                rightDiagUpper[1] = (*v)[vi.index()-2];
                rightDiagUpper[2] = (*v)[vi.index()-1];
                rightDiagUpper[3] = (*v)[vi.index()];
                rightDiagUpper[4] = (*v)[vi.index()+1];
                rightDiagUpper[5] = (*v)[vi.index()+2];

                rightDiagUpper[6] = (*u)[ui.index()-3];
                rightDiagUpper[7] = (*u)[ui.index()-2];
                rightDiagUpper[8] = (*u)[ui.index()-1];
                rightDiagUpper[9] = (*u)[ui.index()];
                rightDiagUpper[10] = (*u)[ui.index()+1];
                rightDiagUpper[11] = (*u)[ui.index()+2];
            }
            DenseVector<WorkPrec_> rightDiagLower(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagLower[0] = (*v)[vi.index()];
                rightDiagLower[1] = (*v)[vi.index()+1];
                rightDiagLower[2] = (*v)[vi.index()+2];
                rightDiagLower[3] = (*v)[vi.index()+3];
                rightDiagLower[4] = (*v)[vi.index()+4];
                rightDiagLower[5] = (*v)[vi.index()+5];

                rightDiagLower[6] = (*u)[ui.index()];
                rightDiagLower[7] = (*u)[ui.index()+1];
                rightDiagLower[8] = (*u)[ui.index()+2];
                rightDiagLower[9] = (*u)[ui.index()+3];
                rightDiagLower[10] = (*u)[ui.index()+4];
                rightDiagLower[11] = (*u)[ui.index()+5];
            }

            //Compute ScalarProducts and Theta - value for diagonal (first theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );
            WorkPrec_ thetaXPlus_i_j;
            if(minusLower!=0)
            {
                 thetaXPlus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (first theta):
            WorkPrec_ thetaXPlus_i_j_limited = min_mod_limiter(thetaXPlus_i_j);

            //Second theta (thetaXPlus_iMinus1_j):
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagUpper[0] = (*v)[vi.index()-6];
                rightDiagUpper[1] = (*v)[vi.index()-5];
                rightDiagUpper[2] = (*v)[vi.index()-4];
                rightDiagUpper[3] = (*v)[vi.index()-3];
                rightDiagUpper[4] = (*v)[vi.index()-2];
                rightDiagUpper[5] = (*v)[vi.index()-1];

                rightDiagUpper[6] = (*u)[ui.index()-6];
                rightDiagUpper[7] = (*u)[ui.index()-5];
                rightDiagUpper[8] = (*u)[ui.index()-4];
                rightDiagUpper[9] = (*u)[ui.index()-3];
                rightDiagUpper[10] = (*u)[ui.index()-2];
                rightDiagUpper[11] = (*u)[ui.index()-1];
                



            }
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagLower[0] = (*v)[vi.index()-3];
                rightDiagLower[1] = (*v)[vi.index()-2];
                rightDiagLower[2] = (*v)[vi.index()-1];
                rightDiagLower[3] = (*v)[vi.index()];
                rightDiagLower[4] = (*v)[vi.index()+1];
                rightDiagLower[5] = (*v)[vi.index()+2];

                rightDiagLower[6] = (*u)[ui.index()-3];
                rightDiagLower[7] = (*u)[ui.index()-2];
                rightDiagLower[8] = (*u)[ui.index()-1];
                rightDiagLower[9] = (*u)[ui.index()];
                rightDiagLower[10] = (*u)[ui.index()+1];
                rightDiagLower[11] = (*u)[ui.index()+2];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (second theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );

            if(minusLower!=0)
            {
                 thetaXPlus_iMinus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_iMinus1_j = 0;
            }
            //Compute limitation of Theta for diagonal (second theta):
            thetaXPlus_iMinus1_j_limited = min_mod_limiter(thetaXPlus_iMinus1_j);

            //Third theta (thetaXMinus_i_j):
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagUpper[0] = (*v)[vi.index()-3];
                rightDiagUpper[1] = (*v)[vi.index()-2];
                rightDiagUpper[2] = (*v)[vi.index()-1];
                rightDiagUpper[3] = (*v)[vi.index()];
                rightDiagUpper[4] = (*v)[vi.index()+1];
                rightDiagUpper[5] = (*v)[vi.index()+2];

                rightDiagUpper[6] = (*u)[ui.index()-3];
                rightDiagUpper[7] = (*u)[ui.index()-2];
                rightDiagUpper[8] = (*u)[ui.index()-1];
                rightDiagUpper[9] = (*u)[ui.index()];
                rightDiagUpper[10] = (*u)[ui.index()+1];
                rightDiagUpper[11] = (*u)[ui.index()+2];
                



            }
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightDiagLower[0] = (*v)[vi.index()];
                rightDiagLower[1] = (*v)[vi.index()+1];
                rightDiagLower[2] = (*v)[vi.index()+2];
                rightDiagLower[3] = (*v)[vi.index()+3];
                rightDiagLower[4] = (*v)[vi.index()+4];
                rightDiagLower[5] = (*v)[vi.index()+5];

                rightDiagLower[6] = (*u)[ui.index()];
                rightDiagLower[7] = (*u)[ui.index()+1];
                rightDiagLower[8] = (*u)[ui.index()+2];
                rightDiagLower[9] = (*u)[ui.index()+3];
                rightDiagLower[10] = (*u)[ui.index()+4];
                rightDiagLower[11] = (*u)[ui.index()+5];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (third theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );
            WorkPrec_ thetaXMinus_i_j;
            if(minusLower!=0)
            {
                 thetaXMinus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (third theta):
            WorkPrec_ thetaXMinus_i_j_limited = min_mod_limiter(thetaXMinus_i_j);

            //Compute matrix element value:
            m1diag[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXPlus_i_j_limited + thetaXPlus_iMinus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            m3diag[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXPlus_i_j_limited + thetaXPlus_iMinus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            //FINISHED diagonal.

            //Prepare right operands for band_+1:
            //First Theta - value (thetaXPlus_i_j)
            DenseVector<WorkPrec_> rightPlus1Upper(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Upper[0] = (*v)[vi.index()-3];
                rightPlus1Upper[1] = (*v)[vi.index()-2];
                rightPlus1Upper[2] = (*v)[vi.index()-1];
                rightPlus1Upper[3] = (*v)[vi.index()];
                rightPlus1Upper[4] = (*v)[vi.index()+1];
                rightPlus1Upper[5] = (*v)[vi.index()+2];

                rightPlus1Upper[6] = (*u)[ui.index()-3];
                rightPlus1Upper[7] = (*u)[ui.index()-2];
                rightPlus1Upper[8] = (*u)[ui.index()-1];
                rightPlus1Upper[9] = (*u)[ui.index()];
                rightPlus1Upper[10] = (*u)[ui.index()+1];
                rightPlus1Upper[11] = (*u)[ui.index()+2];
            }
            DenseVector<WorkPrec_> rightPlus1Lower(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Lower[0] = (*v)[vi.index()];
                rightPlus1Lower[1] = (*v)[vi.index()+1];
                rightPlus1Lower[2] = (*v)[vi.index()+2];
                rightPlus1Lower[3] = (*v)[vi.index()+3];
                rightPlus1Lower[4] = (*v)[vi.index()+4];
                rightPlus1Lower[5] = (*v)[vi.index()+5];

                rightPlus1Lower[6] = (*u)[ui.index()];
                rightPlus1Lower[7] = (*u)[ui.index()+1];
                rightPlus1Lower[8] = (*u)[ui.index()+2];
                rightPlus1Lower[9] = (*u)[ui.index()+3];
                rightPlus1Lower[10] = (*u)[ui.index()+4];
                rightPlus1Lower[11] = (*u)[ui.index()+5];
            }

            //Compute ScalarProducts and Theta - value for diagonal (first theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );
            thetaXPlus_i_j;
            if(minusLower!=0)
            {
                 thetaXPlus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (first theta):
            thetaXPlus_i_j_limited = min_mod_limiter(thetaXPlus_i_j);

            //Second theta (thetaXMinus_iPlus1_j):
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Upper[0] = (*v)[vi.index()];
                rightPlus1Upper[1] = (*v)[vi.index()+1];
                rightPlus1Upper[2] = (*v)[vi.index()+2];
                rightPlus1Upper[3] = (*v)[vi.index()+3];
                rightPlus1Upper[4] = (*v)[vi.index()+4];
                rightPlus1Upper[5] = (*v)[vi.index()+5];

                rightPlus1Upper[6] = (*u)[ui.index()];
                rightPlus1Upper[7] = (*u)[ui.index()+1];
                rightPlus1Upper[8] = (*u)[ui.index()+2];
                rightPlus1Upper[9] = (*u)[ui.index()+3];
                rightPlus1Upper[10] = (*u)[ui.index()+4];
                rightPlus1Upper[11] = (*u)[ui.index()+5];
                



            }
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Lower[0] = (*v)[vi.index()+3];
                rightPlus1Lower[1] = (*v)[vi.index()+4];
                rightPlus1Lower[2] = (*v)[vi.index()+5];
                rightPlus1Lower[3] = (*v)[vi.index()+6];
                rightPlus1Lower[4] = (*v)[vi.index()+7];
                rightPlus1Lower[5] = (*v)[vi.index()+8];

                rightPlus1Lower[6] = (*u)[ui.index()+3];
                rightPlus1Lower[7] = (*u)[ui.index()+4];
                rightPlus1Lower[8] = (*u)[ui.index()+5];
                rightPlus1Lower[9] = (*u)[ui.index()+6];
                rightPlus1Lower[10] = (*u)[ui.index()+7];
                rightPlus1Lower[11] = (*u)[ui.index()+8];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (second theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );
            WorkPrec_ thetaXMinus_iPlus1_j;
            if(minusLower!=0)
            {
                 thetaXMinus_iPlus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_iPlus1_j = 0;
            }
            //Compute limitation of Theta for diagonal (second theta):
            WorkPrec_ thetaXMinus_iPlus1_j_limited = min_mod_limiter(thetaXMinus_iPlus1_j);

            //Third theta (thetaXMinus_i_j):
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Upper[0] = (*v)[vi.index()-3];
                rightPlus1Upper[1] = (*v)[vi.index()-2];
                rightPlus1Upper[2] = (*v)[vi.index()-1];
                rightPlus1Upper[3] = (*v)[vi.index()];
                rightPlus1Upper[4] = (*v)[vi.index()+1];
                rightPlus1Upper[5] = (*v)[vi.index()+2];

                rightPlus1Upper[6] = (*u)[ui.index()-3];
                rightPlus1Upper[7] = (*u)[ui.index()-2];
                rightPlus1Upper[8] = (*u)[ui.index()-1];
                rightPlus1Upper[9] = (*u)[ui.index()];
                rightPlus1Upper[10] = (*u)[ui.index()+1];
                rightPlus1Upper[11] = (*u)[ui.index()+2];
                



            }
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus1Lower[0] = (*v)[vi.index()];
                rightPlus1Lower[1] = (*v)[vi.index()+1];
                rightPlus1Lower[2] = (*v)[vi.index()+2];
                rightPlus1Lower[3] = (*v)[vi.index()+3];
                rightPlus1Lower[4] = (*v)[vi.index()+4];
                rightPlus1Lower[5] = (*v)[vi.index()+5];

                rightPlus1Lower[6] = (*u)[ui.index()];
                rightPlus1Lower[7] = (*u)[ui.index()+1];
                rightPlus1Lower[8] = (*u)[ui.index()+2];
                rightPlus1Lower[9] = (*u)[ui.index()+3];
                rightPlus1Lower[10] = (*u)[ui.index()+4];
                rightPlus1Lower[11] = (*u)[ui.index()+5];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (third theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );

            if(minusLower!=0)
            {
                 thetaXMinus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (third theta):
            thetaXMinus_i_j_limited = min_mod_limiter(thetaXMinus_i_j);

            //Compute matrix element value:
            m1bandPlus1[ui.index()] = -(_delta_t / 4*_delta_x) * (2 + thetaXPlus_i_j_limited + thetaXMinus_iPlus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            m3bandPlus1[ui.index()] = -(_delta_t / 4*_delta_x) * (2 + thetaXPlus_i_j_limited + thetaXMinus_iPlus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            //FINISHED band_+1.

            //band_+2 (thetaXMinus_iPlus1_j):
            DenseVector<WorkPrec_> rightPlus2Upper(ulint(12), ulint(0));
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus2Upper[0] = (*v)[vi.index()];
                rightPlus2Upper[1] = (*v)[vi.index()+1];
                rightPlus2Upper[2] = (*v)[vi.index()+2];
                rightPlus2Upper[3] = (*v)[vi.index()+3];
                rightPlus2Upper[4] = (*v)[vi.index()+4];
                rightPlus2Upper[5] = (*v)[vi.index()+5];

                rightPlus2Upper[6] = (*u)[ui.index()];
                rightPlus2Upper[7] = (*u)[ui.index()+1];
                rightPlus2Upper[8] = (*u)[ui.index()+2];
                rightPlus2Upper[9] = (*u)[ui.index()+3];
                rightPlus2Upper[10] = (*u)[ui.index()+4];
                rightPlus2Upper[11] = (*u)[ui.index()+5];
                



            }
            DenseVector<WorkPrec_> rightPlus2Lower(ulint(12), ulint(0));
 
            if(ui.index()>=6 && ui.index()<=ui_END.index()-6)
            {
                rightPlus2Lower[0] = (*v)[vi.index()+3];
                rightPlus2Lower[1] = (*v)[vi.index()+4];
                rightPlus2Lower[2] = (*v)[vi.index()+5];
                rightPlus2Lower[3] = (*v)[vi.index()+6];
                rightPlus2Lower[4] = (*v)[vi.index()+7];
                rightPlus2Lower[5] = (*v)[vi.index()+8];

                rightPlus2Lower[6] = (*u)[ui.index()+3];
                rightPlus2Lower[7] = (*u)[ui.index()+4];
                rightPlus2Lower[8] = (*u)[ui.index()+5];
                rightPlus2Lower[9] = (*u)[ui.index()+6];
                rightPlus2Lower[10] = (*u)[ui.index()+7];
                rightPlus2Lower[11] = (*u)[ui.index()+8];
                
            }
            //Compute ScalarProducts and Theta - value for band_+2:
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );

            if(minusLower!=0)
            {
                 thetaXMinus_iPlus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_iPlus1_j = 0;
            }
            //Compute limitation of Theta for band_+2:
            thetaXMinus_iPlus1_j_limited = min_mod_limiter(thetaXMinus_iPlus1_j);

            //Compute matrix element value:
            m1bandPlus2[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXMinus_iPlus1_j_limited);
            m3bandPlus2[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXMinus_iPlus1_j_limited);
            //FINISHED band_+2.

            //Iterate:
            ++ui;++vi;
        }
        m1.insert_band(0, m1diag.copy());
        m1.insert_band(3, m1bandPlus1.copy());
        m1.insert_band(6, m1bandPlus2.copy());
        m1.insert_band(-3,m1bandMinus1.copy());
        m3.insert_band(0, m3diag.copy());
        m3.insert_band(3, m3bandPlus1.copy());
        m3.insert_band(6, m3bandPlus2.copy());
        m3.insert_band(-3,m3bandMinus1.copy());
 

        cout << "M_1:" << stringify(m1.band(ulint(0))) << endl;
        cout << "M_1:" << stringify(m1.band(ulint(3))) << endl;
        cout << "M_1:" << stringify(m1.band(ulint(6))) << endl;
        cout << "M_1:" << stringify(m1.band(ulint(-3))) << endl;
        std::cout << "Finished Matrix Assembly 1.\n";
 
    }
    /**
     *
     * Matrix assembly. DEBUG version.
     *
     * \param  m2 The first matrix to assemble.
     * \param  m4 The second matrix to assemble.
     *
     **/
    template<typename ResPrec_ , typename PredictionPrec1_, typename PredictionPrec2_, typename InitPrec1_, typename InitPrec2_>
    template<typename WorkPrec_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>
        ::_assemble_matrix2_DEBUG(BandedMatrix<WorkPrec_>& m2, BandedMatrix<WorkPrec_>& m4, DenseVector<WorkPrec_>* u, DenseVector<WorkPrec_>* v)
    {
        ///The bands containing data.
        /*DenseVector<WorkPrec_> m2diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m2bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m2bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m2bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m4diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m4bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m4bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m4bandMinus1(_u->size(),ulint( 0));//three
        m2.band(ulint(0)) = m2diag;
        m2.band(ulint(3*(_d_width-1))) = m2bandPlus1;
        m2.band(ulint(6*(_d_width-1))) = m2bandPlus2;
        m2.band(ulint(-3*(_d_width-1))) = m2bandMinus1;
        m4.band(ulint(0)) = m4diag;
        m4.band(ulint(3*(_d_width-1))) = m4bandPlus1;
        m4.band(ulint(6*(_d_width-1))) = m4bandPlus2;
        m4.band(ulint(-3*(_d_width-1))) = m4bandMinus1;
        */
        DenseVector<WorkPrec_> m2diag(_u->size(), ulint(0));      //zero
        DenseVector<WorkPrec_> m2bandPlus1(_u->size(), ulint(0)); //one
        DenseVector<WorkPrec_> m2bandPlus2(_u->size(), ulint(0)); //two
        DenseVector<WorkPrec_> m2bandMinus1(_u->size(), ulint(0));//three
        DenseVector<WorkPrec_> m4diag(_u->size(),ulint( 0));      //zero
        DenseVector<WorkPrec_> m4bandPlus1(_u->size(),ulint (0)); //one
        DenseVector<WorkPrec_> m4bandPlus2(_u->size(),ulint (0)); //two
        DenseVector<WorkPrec_> m4bandMinus1(_u->size(),ulint( 0));//three

        typename DenseVector<WorkPrec_>::ElementIterator ui(u->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator ui_END(u->end_elements());
        typename DenseVector<WorkPrec_>::ElementIterator vi(u->begin_elements());
        typename DenseVector<WorkPrec_>::ElementIterator vi_END(u->end_elements());
        
        //Preparing left operand for all scalar products
        DenseVector<WorkPrec_> constantVector(ulint(12), ulint(0));
        constantVector[0] = WorkPrec_(-1);
        constantVector[1] = WorkPrec_(-1);
        constantVector[2] = WorkPrec_(-1);
        constantVector[3] = WorkPrec_(1);
        constantVector[4] = WorkPrec_(1);
        constantVector[5] = WorkPrec_(1);
        constantVector[6] = WorkPrec_(-((*_c)[0]));
        constantVector[7] = WorkPrec_(-((*_c)[1]));
        constantVector[8] = WorkPrec_(-((*_c)[2]));
        constantVector[9]  = WorkPrec_((*_c)[0]);
        constantVector[10] = WorkPrec_((*_c)[1]);
        constantVector[11] = WorkPrec_((*_c)[2]);
        //Preparing left operand for all scalar products
        DenseVector<WorkPrec_> constantVectorMinus(ulint(12), ulint(0));
        constantVectorMinus[0] = WorkPrec_(-1);
        constantVectorMinus[1] = WorkPrec_(-1);
        constantVectorMinus[2] = WorkPrec_(-1);
        constantVectorMinus[3] = WorkPrec_(1);
        constantVectorMinus[4] = WorkPrec_(1);
        constantVectorMinus[5] = WorkPrec_(1);
        constantVectorMinus[6] = WorkPrec_((*_c)[0]);
        constantVectorMinus[7] = WorkPrec_((*_c)[1]);
        constantVectorMinus[8] = WorkPrec_((*_c)[2]);
        constantVectorMinus[9]  = WorkPrec_(-(*_c)[0]);
        constantVectorMinus[10] = WorkPrec_(-(*_c)[1]);
        constantVectorMinus[11] = WorkPrec_(-(*_c)[2]);

 
        while(ui!=ui_END)
        {
            //Prepare right operands for band_-1:
            //thetaXPlus_iMinus1_j
            DenseVector<WorkPrec_> rightMinus1Upper(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightMinus1Upper[0] = (*v)[vi.index()-6*(_d_width+4)];
                rightMinus1Upper[1] = (*v)[vi.index()-6*(_d_width+4)+1];
                rightMinus1Upper[2] = (*v)[vi.index()-6*(_d_width+4)+2];
                rightMinus1Upper[3] = (*v)[vi.index()-3*(_d_width+4)];
                rightMinus1Upper[4] = (*v)[vi.index()-3*(_d_width+4)+1];
                rightMinus1Upper[5] = (*v)[vi.index()-3*(_d_width+4)+2];

                rightMinus1Upper[6] = (*u)[ui.index()-6*(_d_width+4)];
                rightMinus1Upper[7] = (*u)[ui.index()-6*(_d_width+4)+1];
                rightMinus1Upper[8] = (*u)[ui.index()-6*(_d_width+4)+2];
                rightMinus1Upper[9] = (*u)[ui.index()-3*(_d_width+4)];
                rightMinus1Upper[10] = (*u)[ui.index()-3*(_d_width+4)+1];
                rightMinus1Upper[11] = (*u)[ui.index()-3*(_d_width+4)+2];
                



            }
            DenseVector<WorkPrec_> rightMinus1Lower(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4)) 
            {
                rightMinus1Lower[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightMinus1Lower[1] = (*v)[vi.index()-3*(_d_width+4)+1];
                rightMinus1Lower[2] = (*v)[vi.index()-3*(_d_width+4)+2];
                rightMinus1Lower[3] = (*v)[vi.index()];
                rightMinus1Lower[4] = (*v)[vi.index()+1];
                rightMinus1Lower[5] = (*v)[vi.index()+2];

                rightMinus1Lower[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightMinus1Lower[7] = (*u)[ui.index()-3*(_d_width+4)+1];
                rightMinus1Lower[8] = (*u)[ui.index()-3*(_d_width+4)+2];
                rightMinus1Lower[9] = (*u)[ui.index()];
                rightMinus1Lower[10] = (*u)[ui.index()+1];
                rightMinus1Lower[11] = (*u)[ui.index()+2];
                
            }

            //Compute ScalarProducts and Theta - value for band_-1:
            WorkPrec_ minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightMinus1Upper, constantVector );
            WorkPrec_ minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightMinus1Lower, constantVector );
            WorkPrec_ thetaXPlus_iMinus1_j;
            if(minusLower!=0)
            {
                 thetaXPlus_iMinus1_j = minusUpper / minusLower;
            }
            //Compute limitation of Theta for band_-1:
            WorkPrec_ thetaXPlus_iMinus1_j_limited = min_mod_limiter(thetaXPlus_iMinus1_j);

            //Compute matrix element value:
            m2bandMinus1[ui.index()] = (_delta_t / 4*_delta_x) * (2 -  thetaXPlus_iMinus1_j_limited);
            m4bandMinus1[ui.index()] = (_delta_t / 4*_delta_x) * (2 -  thetaXPlus_iMinus1_j_limited);
            //FINISHED band_-1.

            //Prepare right operands for diagonal:
            //First Theta - value (thetaXPlus_i_j)
            DenseVector<WorkPrec_> rightDiagUpper(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagUpper[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagUpper[1] = (*v)[vi.index()-3*(_d_width+4)+1];
                rightDiagUpper[2] = (*v)[vi.index()-3*(_d_width+4)+2];
                rightDiagUpper[3] = (*v)[vi.index()];
                rightDiagUpper[4] = (*v)[vi.index()+1];
                rightDiagUpper[5] = (*v)[vi.index()+2];

                rightDiagUpper[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagUpper[7] = (*u)[ui.index()-3*(_d_width+4)+1];
                rightDiagUpper[8] = (*u)[ui.index()-3*(_d_width+4)+2];
                rightDiagUpper[9] = (*u)[ui.index()];
                rightDiagUpper[10] = (*u)[ui.index()+1];
                rightDiagUpper[11] = (*u)[ui.index()+2];
            }
            DenseVector<WorkPrec_> rightDiagLower(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagLower[0] = (*v)[vi.index()];
                rightDiagLower[1] = (*v)[vi.index()+1];
                rightDiagLower[2] = (*v)[vi.index()+2];
                rightDiagLower[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightDiagLower[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightDiagLower[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightDiagLower[6] = (*u)[ui.index()];
                rightDiagLower[7] = (*u)[ui.index()+1];
                rightDiagLower[8] = (*u)[ui.index()+2];
                rightDiagLower[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightDiagLower[10] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightDiagLower[11] = (*u)[ui.index()+3*(_d_width+4)+2];
            }

            //Compute ScalarProducts and Theta - value for diagonal (first theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );
            WorkPrec_ thetaXPlus_i_j;
            if(minusLower!=0)
            {
                 thetaXPlus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (first theta):
            WorkPrec_ thetaXPlus_i_j_limited = min_mod_limiter(thetaXPlus_i_j);

            //Second theta (thetaXPlus_iMinus1_j):
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagUpper[0] = (*v)[vi.index()-6*(_d_width)];
                rightDiagUpper[1] = (*v)[vi.index()-6*(_d_width)+1];
                rightDiagUpper[2] = (*v)[vi.index()-6*(_d_width)+2];
                rightDiagUpper[3] = (*v)[vi.index()-3*(_d_width)];
                rightDiagUpper[4] = (*v)[vi.index()-3*(_d_width)+1];
                rightDiagUpper[5] = (*v)[vi.index()-3*(_d_width)+2];

                rightDiagUpper[6] = (*u)[ui.index()-6*(_d_width)];
                rightDiagUpper[7] = (*u)[ui.index()-6*(_d_width)+1];
                rightDiagUpper[8] = (*u)[ui.index()-6*(_d_width)+2];
                rightDiagUpper[9] = (*u)[ui.index()-3*_d_width];
                rightDiagUpper[10] = (*u)[ui.index()-3*(_d_width)+1];
                rightDiagUpper[11] = (*u)[ui.index()-3*(_d_width)+2];
                



            }
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagLower[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagLower[1] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagLower[2] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagLower[3] = (*v)[vi.index()];
                rightDiagLower[4] = (*v)[vi.index()+1];
                rightDiagLower[5] = (*v)[vi.index()+2];

                rightDiagLower[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagLower[7] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagLower[8] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagLower[9] = (*u)[ui.index()];
                rightDiagLower[10] = (*u)[ui.index()+1];
                rightDiagLower[11] = (*u)[ui.index()+2];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (second theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );

            if(minusLower!=0)
            {
                 thetaXPlus_iMinus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_iMinus1_j = 0;
            }
            //Compute limitation of Theta for diagonal (second theta):
            thetaXPlus_iMinus1_j_limited = min_mod_limiter(thetaXPlus_iMinus1_j);

            //Third theta (thetaXMinus_i_j):
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagUpper[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagUpper[1] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagUpper[2] = (*v)[vi.index()-3*(_d_width+4)];
                rightDiagUpper[3] = (*v)[vi.index()];
                rightDiagUpper[4] = (*v)[vi.index()+1];
                rightDiagUpper[5] = (*v)[vi.index()+2];

                rightDiagUpper[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagUpper[7] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagUpper[8] = (*u)[ui.index()-3*(_d_width+4)];
                rightDiagUpper[9] = (*u)[ui.index()];
                rightDiagUpper[10] = (*u)[ui.index()+1];
                rightDiagUpper[11] = (*u)[ui.index()+2];
                



            }
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightDiagLower[0] = (*v)[vi.index()];
                rightDiagLower[1] = (*v)[vi.index()+1];
                rightDiagLower[2] = (*v)[vi.index()+2];
                rightDiagLower[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightDiagLower[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightDiagLower[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightDiagLower[6] = (*u)[ui.index()];
                rightDiagLower[7] = (*u)[ui.index()+1];
                rightDiagLower[8] = (*u)[ui.index()+2];
                rightDiagLower[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightDiagLower[10] = (*u)[ui.index()+3*(_d_width+4)];
                rightDiagLower[11] = (*u)[ui.index()+3*(_d_width+4)];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (third theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );
            WorkPrec_ thetaXMinus_i_j;
            if(minusLower!=0)
            {
                 thetaXMinus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (third theta):
            WorkPrec_ thetaXMinus_i_j_limited = min_mod_limiter(thetaXMinus_i_j);

            //Compute matrix element value:
            m2diag[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXPlus_i_j_limited + thetaXPlus_iMinus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            m4diag[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXPlus_i_j_limited + thetaXPlus_iMinus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            //FINISHED diagonal.

            //Prepare right operands for band_+1:
            //First Theta - value (thetaXPlus_i_j)
            DenseVector<WorkPrec_> rightPlus1Upper(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Upper[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightPlus1Upper[1] = (*v)[vi.index()-3*(_d_width+4)+1];
                rightPlus1Upper[2] = (*v)[vi.index()-3*(_d_width+4)+2];
                rightPlus1Upper[3] = (*v)[vi.index()];
                rightPlus1Upper[4] = (*v)[vi.index()+1];
                rightPlus1Upper[5] = (*v)[vi.index()+2];

                rightPlus1Upper[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightPlus1Upper[7] = (*u)[ui.index()-3*(_d_width+4)+1];
                rightPlus1Upper[8] = (*u)[ui.index()-3*(_d_width+4)+2];
                rightPlus1Upper[9] = (*u)[ui.index()];
                rightPlus1Upper[10] = (*u)[ui.index()+1];
                rightPlus1Upper[11] = (*u)[ui.index()+2];
            }
            DenseVector<WorkPrec_> rightPlus1Lower(ulint(12), ulint(0));
            if(ui.index()>=6 *(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Lower[0] = (*v)[vi.index()];
                rightPlus1Lower[1] = (*v)[vi.index()+1];
                rightPlus1Lower[2] = (*v)[vi.index()+2];
                rightPlus1Lower[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus1Lower[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus1Lower[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightPlus1Lower[6] = (*u)[ui.index()];
                rightPlus1Lower[7] = (*u)[ui.index()+1];
                rightPlus1Lower[8] = (*u)[ui.index()+2];
                rightPlus1Lower[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus1Lower[10] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus1Lower[11] = (*u)[ui.index()+3*(_d_width+4)+2];
            }

            //Compute ScalarProducts and Theta - value for diagonal (first theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );
            thetaXPlus_i_j;
            if(minusLower!=0)
            {
                 thetaXPlus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXPlus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (first theta):
            thetaXPlus_i_j_limited = min_mod_limiter(thetaXPlus_i_j);

            //Second theta (thetaXMinus_iPlus1_j):
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Upper[0] = (*v)[vi.index()];
                rightPlus1Upper[1] = (*v)[vi.index()+1];
                rightPlus1Upper[2] = (*v)[vi.index()+2];
                rightPlus1Upper[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus1Upper[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus1Upper[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightPlus1Upper[6] = (*u)[ui.index()];
                rightPlus1Upper[7] = (*u)[ui.index()+1];
                rightPlus1Upper[8] = (*u)[ui.index()+2];
                rightPlus1Upper[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus1Upper[10] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus1Upper[11] = (*u)[ui.index()+3*(_d_width+4)+2];
                



            }
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Lower[0] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus1Lower[1] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus1Lower[2] = (*v)[vi.index()+3*(_d_width+4)+2];
                rightPlus1Lower[3] = (*v)[vi.index()+6*(_d_width+4)];
                rightPlus1Lower[4] = (*v)[vi.index()+6*(_d_width+4)+1];
                rightPlus1Lower[5] = (*v)[vi.index()+6*(_d_width+4)+2];

                rightPlus1Lower[6] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus1Lower[7] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus1Lower[8] = (*u)[ui.index()+3*(_d_width+4)+2];
                rightPlus1Lower[9] = (*u)[ui.index()+6*(_d_width+4)];
                rightPlus1Lower[10] = (*u)[ui.index()+6*(_d_width+4)+1];
                rightPlus1Lower[11] = (*u)[ui.index()+6*(_d_width+4)+2];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (second theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );
            WorkPrec_ thetaXMinus_iPlus1_j;
            if(minusLower!=0)
            {
                 thetaXMinus_iPlus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_iPlus1_j = 0;
            }
            //Compute limitation of Theta for diagonal (second theta):
            WorkPrec_ thetaXMinus_iPlus1_j_limited = min_mod_limiter(thetaXMinus_iPlus1_j);

            //Third theta (thetaXMinus_i_j):
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Upper[0] = (*v)[vi.index()-3*(_d_width+4)];
                rightPlus1Upper[1] = (*v)[vi.index()-3*(_d_width+4)+1];
                rightPlus1Upper[2] = (*v)[vi.index()-3*(_d_width+4)+2];
                rightPlus1Upper[3] = (*v)[vi.index()];
                rightPlus1Upper[4] = (*v)[vi.index()+1];
                rightPlus1Upper[5] = (*v)[vi.index()+2];

                rightPlus1Upper[6] = (*u)[ui.index()-3*(_d_width+4)];
                rightPlus1Upper[7] = (*u)[ui.index()-3*(_d_width+4)+1];
                rightPlus1Upper[8] = (*u)[ui.index()-3*(_d_width+4)+2];
                rightPlus1Upper[9] = (*u)[ui.index()];
                rightPlus1Upper[10] = (*u)[ui.index()+1];
                rightPlus1Upper[11] = (*u)[ui.index()+2];
                



            }
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus1Lower[0] = (*v)[vi.index()];
                rightPlus1Lower[1] = (*v)[vi.index()+1];
                rightPlus1Lower[2] = (*v)[vi.index()+2];
                rightPlus1Lower[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus1Lower[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus1Lower[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightPlus1Lower[6] = (*u)[ui.index()];
                rightPlus1Lower[7] = (*u)[ui.index()+1];
                rightPlus1Lower[8] = (*u)[ui.index()+2];
                rightPlus1Lower[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus1Lower[10] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus1Lower[11] = (*u)[ui.index()+3*(_d_width+4)+2];
                
            }
            //Compute ScalarProducts and Theta - value for diagonal (third theta):
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVector );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVector );

            if(minusLower!=0)
            {
                 thetaXMinus_i_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_i_j = 0;
            }
            //Compute limitation of Theta for diagonal (third theta):
            thetaXMinus_i_j_limited = min_mod_limiter(thetaXMinus_i_j);

            //Compute matrix element value:
            m2bandPlus1[ui.index()] = -(_delta_t / 4*_delta_x) * (2 + thetaXPlus_i_j_limited + thetaXMinus_iPlus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            m4bandPlus1[ui.index()] = -(_delta_t / 4*_delta_x) * (2 + thetaXPlus_i_j_limited + thetaXMinus_iPlus1_j_limited + 
                                                       thetaXMinus_i_j_limited );
            //FINISHED band_+1.

            //band_+2 (thetaXMinus_iPlus1_j):
            DenseVector<WorkPrec_> rightPlus2Upper(ulint(12), ulint(0));
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus2Upper[0] = (*v)[vi.index()];
                rightPlus2Upper[1] = (*v)[vi.index()+1];
                rightPlus2Upper[2] = (*v)[vi.index()+2];
                rightPlus2Upper[3] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus2Upper[4] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus2Upper[5] = (*v)[vi.index()+3*(_d_width+4)+2];

                rightPlus2Upper[6] = (*u)[ui.index()];
                rightPlus2Upper[7] = (*u)[ui.index()+1];
                rightPlus2Upper[8] = (*u)[ui.index()+2];
                rightPlus2Upper[9] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus2Upper[10] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus2Upper[11] = (*u)[ui.index()+3*(_d_width+4)+2];
                



            }
            DenseVector<WorkPrec_> rightPlus2Lower(ulint(12), ulint(0));
 
            if(ui.index()>=6*(_d_width+4) && ui.index()<=ui_END.index()-6*(_d_width+4))
            {
                rightPlus2Lower[0] = (*v)[vi.index()+3*(_d_width+4)];
                rightPlus2Lower[1] = (*v)[vi.index()+3*(_d_width+4)+1];
                rightPlus2Lower[2] = (*v)[vi.index()+3*(_d_width+4)+2];
                rightPlus2Lower[3] = (*v)[vi.index()+6*(_d_width+4)];
                rightPlus2Lower[4] = (*v)[vi.index()+6*(_d_width+4)+1];
                rightPlus2Lower[5] = (*v)[vi.index()+6*(_d_width+4)+2];

                rightPlus2Lower[6] = (*u)[ui.index()+3*(_d_width+4)];
                rightPlus2Lower[7] = (*u)[ui.index()+3*(_d_width+4)+1];
                rightPlus2Lower[8] = (*u)[ui.index()+3*(_d_width+4)+2];
                rightPlus2Lower[9] = (*u)[ui.index()+6*(_d_width+4)];
                rightPlus2Lower[10] = (*u)[ui.index()+6*(_d_width+4)+1];
                rightPlus2Lower[11] = (*u)[ui.index()+6*(_d_width+4)+2];
                
            }
            //Compute ScalarProducts and Theta - value for band_+2:
            minusUpper = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagUpper, constantVectorMinus );
            minusLower = ScalarProduct<>::value<WorkPrec_, WorkPrec_>( rightDiagLower, constantVectorMinus );

            if(minusLower!=0)
            {
                 thetaXMinus_iPlus1_j = minusUpper / minusLower;
            }
            else
            {
                thetaXMinus_iPlus1_j = 0;
            }
            //Compute limitation of Theta for band_+2:
            thetaXMinus_iPlus1_j_limited = min_mod_limiter(thetaXMinus_iPlus1_j);

            //Compute matrix element value:
            m2bandPlus2[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXMinus_iPlus1_j_limited);
            m4bandPlus2[ui.index()] = (_delta_t / 4*_delta_x) * (thetaXMinus_iPlus1_j_limited);
            //FINISHED band_+2.

            //Iterate:
            ++ui;++vi;
        }
        m2.insert_band(0, m2diag.copy());
        m2.insert_band(3*(_d_width +4), m2bandPlus1.copy());
        m2.insert_band(6*(_d_width +4), m2bandPlus2.copy());
        m2.insert_band(-3*(_d_width +4), m2bandMinus1.copy());
        m4.insert_band(0, m4diag.copy());
        m4.insert_band(3*(_d_width +4), m4bandPlus1.copy());
        m4.insert_band(6*(_d_width +4), m4bandPlus2.copy());
        m4.insert_band(-3*(_d_width +4), m4bandMinus1.copy());
        
        cout << "M_2:" << stringify(m2.band(ulint(0))) << endl;
        cout << "M_2:" << stringify(m2.band(ulint(3*(_d_width +4)))) << endl;
        cout << "M_2:" << stringify(m2.band(ulint(6*(_d_width +4)))) << endl;
        cout << "M_2:" << stringify(m2.band(ulint(-3*(_d_width +4)))) << endl;
        std::cout << "Finished Matrix Assembly 2.\n";
 
    }

    /** Implementation of a simple postprocessing method generating a file for GNUPLOT>splot
     *  to display the heightfield in one timestep. The fileformat is "default", which means, that
     *  for each row, a block with tripels is written to the file: <x, y, value>. The row-blocks 
     *  are separated by empty records within the file. When reading the file without plotting, be careful to
     *  have in mind, that a step of 1 from one row to the other means a step of one in x -direction so we
     *  have the "scalarfield" -point of view and not the "matrix" -point of view.
     *  To simply display the file type
     *      gnuplot
     *      splot "filename"
     *  For more details see
     *      gnuplot
     *      help splot
     *  The filenames are consisting of the string "rsolverfield" and the number of the related timestep. The filesuffix is "dat".
     *
     * \param every In order to not have an overload of files, we do not need, one can specify the number by which _solve_time has to be able to be divided by restlessly - obviously this will slow down the computation, so the common use case will be to set every to _solve_time, so that we render only the data for the last timestep.
     *
     * 
     **/ 
    template<typename ResPrec_ , typename PredictionPrec1_, typename PredictionPrec2_, typename InitPrec1_, typename InitPrec2_>
    void RelaxSolver<ResPrec_, PredictionPrec1_, PredictionPrec2_, InitPrec1_, InitPrec2_>::_do_postprocessing(int every)
    {
        if(_solve_time % every == 0 || _solve_time == 0)
        {
            string filename;
            ofstream file;
            filename = "rsolverfield" + stringify(_solve_time) + ".dat";
            
            file.open(filename.c_str());

            ///For every column, traverse elements and write record to file:
            for(ulint x = 0; x < _d_width; ++x)
            {
                for(ulint y = 0; y < _d_height; ++y)
                {
                    string record = stringify(x) + " " + stringify(y) + " " + stringify((*_height)[y][x]) + "\n";
                    file << record;
                }
                //Create empty record after each row:
                file << "\n";
            }
            file.close();
        }
        cout <<"Finished postprocessing." << endl;
    }
 

}
#endif
