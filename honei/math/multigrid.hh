/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBMATH_GUARD_MULTIGRID_HH
#define LIBMATH_GUARD_MULTIGRID_HH 1

#include<honei/la/dense_vector.hh>
#include<honei/la/banded_matrix_q1.hh>
#include<honei/math/methods.hh>
#include<honei/math/conjugate_gradients.hh>
#include<honei/math/jacobi.hh>
#include<honei/la/norm.hh>
#include<honei/la/sum.hh>
#include<honei/la/difference.hh>
#include<honei/la/scaled_sum.hh>
#include<honei/la/scale.hh>
#include<honei/la/product.hh>
#include<honei/la/algorithm.hh>
#include<honei/math/restriction.hh>
#include<honei/math/prolongation.hh>
#include<honei/math/endian_swap.hh>
#include<honei/math/defect.hh>
#include<vector>
#include<string>
#include<fstream>
#include<honei/util/time_stamp.hh>
#include <cmath>

/**
 * Use SOLVER_VERBOSE for detailed iteration-wise output of defect norms.
 */
//#define SOLVER_VERBOSE 1
//#define SOLVER_BENCHMARK


using namespace methods;
namespace honei
{

    /**
     * \brief Configuration and data structure for the Multigrid solver.
     *
     * MGInfo contains all data and configuration parameters needed by the multigrid solvers.
     *
     */
    template <typename Prec_, typename DataLayout_>
    struct MGInfo
        {
            public:
                //configuration constants:

                /**
                 * Defines, wether the solver is used in smoother-mode or not.
                 */
                bool is_smoother;

                /**
                 * Contains constants to encode boundary condition type on the four edges and
                 * the four corners of the parameterplane. Allowed values are 1 vor NEUMANN boundaries
                 * and 2 for DIRICHLET boundaries. Used by Restriction and Prolongation.
                 */
                DenseVector<unsigned long>* macro_border_mask;

                /**
                 * The coasest and finest levels in the cycle.
                 */
                unsigned long min_level;
                unsigned long max_level;

                /**
                 * The maximum number of iterations.
                 */
                unsigned long  n_max_iter;

                /**
                 * Defines, wether the initial guess is the zero-vector or not.
                 */
                bool initial_zero;

                /**
                 * Parameters used for convergence check.
                 */
                Prec_ tolerance;
                bool convergence_check;

                /**
                 * Stores the number of iterations of the MG kernel used by the mixed
                 * precision framework.
                 */
                unsigned long inner_iterations;

                /**
                 * Total number of iterations of the pre- and postsmoother and maximum number
                 * of iterations of the coarse grid solver.
                 */
                unsigned long  n_pre_smooth;
                unsigned long  n_post_smooth;
                unsigned long  n_max_iter_coarse;

                /**
                 * Parameter used for convergence check of the coarse grid solver.
                 */
                Prec_ tolerance_coarse;

                /**
                 * Parameter used for adaptive correction. See kernel code for details.
                 */
                Prec_ adapt_correction_factor;

                /**
                 * STL vectors to store the data on all levels.
                 */
                std::vector<DenseVector<Prec_> > c;
                std::vector<DenseVector<Prec_> > d;
                std::vector<DenseVector<Prec_> > rhs;
                std::vector<DenseVector<Prec_> > x;
                std::vector<DenseVector<Prec_> > diags_inverted;

                std::vector<DataLayout_ > a;
                std::vector<DataLayout_ > prolmats;
                std::vector<DataLayout_ > resmats;
        };

    template<typename Tag_, typename OuterTag_, typename ProlType_, typename SmootherType_, typename CycleType_, typename Mode_>
    struct Multigrid
    {
    };

    template<typename Tag_, typename OuterTag_, typename ProlType_>
    class Multigrid<Tag_, OuterTag_, ProlType_, JAC, CYCLE::V, FIXED>
    {
        template <typename Prec_>
        static DenseVector<Prec_> _multigrid_kernel(BandedMatrixQ1<Prec_>&  system, DenseVector<Prec_>& right_hand_side, unsigned long max_levels, Prec_ * cappa, MGInfo<Prec_, BandedMatrixQ1<Prec_> > & info)
        {
            bool restriction_started(false);
            // compute initial defect
            // when start vector is zero: set d = rhs (we can save one matvec here)
            // when start vector is not zero: d = rhs - A*x
            if(info.initial_zero)
            {
                info.d[info.max_level] = info.rhs[info.max_level];
            }
            else
            {
                DenseVector<Prec_> defect(Defect<Tag_>::value(info.rhs[info.max_level], info.a[info.max_level], info.x[info.max_level]));
                info.d[info.max_level] = defect;
            }

            Prec_ defect, initial_defect;
            // compute norm of initial defect, also used for relative convergence control
            if(!info.is_smoother)
            {
                defect = Norm<vnt_l_two, true, Tag_>::value(info.d[info.max_level]);
                initial_defect = defect;
            }
            else
            {
                defect = Prec_(1e8);
                initial_defect = defect;
            }
            // check if nothing needs to be done
            if(info.convergence_check && defect <= info.tolerance)
            {
                //do nothing
            }
            else
            {
                //start cycles
                unsigned long iter, current_level(info.max_level);
                unsigned long local_cycle[max_levels];
                for(iter = 1 ; iter <= info.n_max_iter ; ++iter)
                {
#ifdef SOLVER_VERBOSE
                    std::cout << iter << "th iteration" <<std::endl;
#endif
                    // set current level to the maximal level
                    current_level = info.max_level;

                    // Initialise the array providing the local cycle information
                    //V cycle!
                    for(unsigned long i(info.min_level + 1); i <= info.max_level; ++i)
                    {
                        local_cycle[i] = 1; //Using 1 to code local V cycle
                    }
                    //cycle loop
                    while(true)
                    {
                        //Restriction loop:
                        // set a flag that the restriction has just started
                        restriction_started = true;
                        while(true)
                        {
                            // for the case that we have actually only one level, do no presmoothing
                            // and jump to the coarse grid solver directly.
                            if (info.min_level == info.max_level)
                                goto endRestrictionLoop;

                            // -------------
                            // presmoothing
                            // -------------

                            if (restriction_started)
                            {
                                // When the restriction loop just started

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], (info.c[current_level]), Prec_(0.7), info.diags_inverted[current_level]);

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], (info.c[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                                Sum<Tag_>::value(info.x[current_level], info.c[current_level]);
                            }
                            else
                            {
                                // otherwise the solution process can be started directly with
                                // the cleared solution vector (as the solution vector itself represents
                                // the defect correction here)
                                //info.x[current_level] = (Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], Prec_(0.7), info.diags_inverted[current_level]);
                                //info.x[current_level] = (Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]),info.x[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                                //END NEWTEST
                            }
                            DenseVector<Prec_> defect_2(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_2;

#ifdef SOLVER_VERBOSE
                            std::cout << "-----------------------------------------------------" << std::endl;
                            std::cout << "Presmoothing ||D|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            //----------------------------------
                            //restriction ("go down one level")
                            //----------------------------------

                            --current_level;

                            // restrict defect from level icurrentLevel+1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the restricted defect vector
                            // depending on Dirichlet mask (see routine for details), and store a copy in RHS

                            /*info.d[current_level] =*/ Restriction<Tag_>::value((info.d[current_level]), (info.d[current_level + 1]), *info.macro_border_mask);
#ifdef SOLVER_VERBOSE
                            std::cout << "Restricted." << std::endl;
#endif
                            info.rhs[current_level] =(info.d[current_level]);

#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            // if we reached the coarsest level exit the restricition loop
                            if (current_level == info.min_level)
                                goto endRestrictionLoop;

                            restriction_started = false;
                        }
endRestrictionLoop:
                        // ----------------------
                        // coarse grid correction
                        // ----------------------

                        if (info.min_level == info.max_level)
                        {
                            // For the case we actually have only one MG level, only
                            // the following coarse grid correction (and no smoothing) is done.

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());

                            DenseVector<Prec_> defect_3(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_3;

                            //
                            // the MG cycle can be exited immediately
                            //
                            goto endCycleLoop;
                        }
                        else
                        {
                            // Otherwise this is a "real" coarse grid correction, which is
                            // started with a zero start vector

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());
#ifdef SOLVER_VERBOSE
                            std::cout << "Coarse Grid solver." << std::endl;
#endif
                        }

                        //-------------
                        //prolongation
                        //-------------
                        while (true)
                        {

                            //------------------------------
                            //prolongation ("go up one level")
                            //------------------------------

                            ++current_level;

                            //
                            // prolongate solution vector from level icurrentLevel-1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the prolongated correction vector
                            // depending on Dirichlet mask passed in from FEAST (see code for details)
                            //
                            /*if (current_level == 4)
                              {
                              info.c[current_level].lock(lm_write_only);
                              info.x[current_level - 1].lock(lm_read_only);
                              info.c[current_level] = Prolongation<tags::CPU>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask);
                              info.c[current_level].unlock(lm_write_only);
                              info.x[current_level - 1].unlock(lm_read_only);
                              }
                              else*/
                            {
                                /*info.c[current_level] =*/ Prolongation<Tag_, ProlType_>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask, info.prolmats[current_level]);
                            }
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongated." << std::endl;
#endif
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||c|| " << Norm<vnt_l_two, true, Tag_>::value(info.c[current_level]) << std::endl;
#endif
                            //
                            // perform adaptive coarse grid correction if required
                            //
                            Prec_ alpha;
                            /*if (info.adapt_correction_factor == 0.0)
                              {
                            //Compute dalpha = (d,c)/(Ac,c) where d is the current residual
                            //(data->d[icurrentLevel]) and c the prolongated correction (c[icurrentLevel])
                            Prec_ d1,d2;
                            d1 = DotProduct<Tag_>::value((info.d[current_level]), (info.c[current_level]));
                            DenseVector<Prec_> t(Product<Tag_>::value((info.a[current_level]), (info.c[current_level])));

                            d2 = DotProduct<Tag_>::value(t, (info.c[current_level]));

                            alpha = Prec_(d1 / d2);
                            }
                            else*/
                            {
                                alpha = info.adapt_correction_factor;
                            }

                            //
                            // add the prolongated correction to the iteration vector
                            // on the current level
                            //
                            ScaledSum<Tag_>::value((info.x[current_level]), (info.c[current_level]), alpha);

#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||x|| " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //--------------
                            //postsmoothing
                            //--------------
                            //
                            // smooth A*x = rhs based on the RHS for that level we stored during restriction

                            //(info.x[current_level]) =(Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]));
                            Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]),(info.x[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]);
                            //end NEWTEST
#ifdef SOLVER_VERBOSE
                            std::cout << "Postsmoothing ||X|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //
                            // update defect
                            //
                            DenseVector<Prec_> defect_4(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_4;
#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
                            std::cout << "-----------------------------------------------------" << std::endl;
#endif
                            // if the maximal level is reached then the MG cycle is finished,
                            // so exit the MG cycle loop
                            if (current_level == info.max_level)
                                goto endCycleLoop;

                            // If the prolongation did not reach the finest level yet, then
                            // determine how to proceed now. Dependent on the type of MG cycle the next
                            // step is a prolongation ("going up") OR a restriction ("going down").
                            if (local_cycle[current_level] == 1)
                            {
                                // local V-cycle --> go on with prolongation ("going up")
                                // this is the case when we are in a global V-cycle or at the end of
                                // a local W-cycle
                                if (true)
                                {
                                    // When we are in a global V- or F-cycle then do a V-cycle on this level
                                    // next time (in case of a global F-cycle (which does a local W-cycle first
                                    // and then a local V-cycle) it means that the local W-cycle
                                    // on this level is finished, such that next time a local V-cycle has to
                                    // be performed)
                                    local_cycle[current_level] = 1;
                                }
                                else
                                {
                                    // In case of a global W-cycle also set the local cycle to a W-cycle
                                    local_cycle[current_level] = 2;
                                }
                            }
                            else
                            {
                                // local W-cycle --> stop prolongation and do restriction again ("going down")
                                // local W-cycle becomes local V-cycle to indicate that half the W-cycle
                                // has been performed
                                local_cycle[current_level] = 1;
                                // exit prolongation loop
                                goto endProlongationLoop;
                            }
                        }
endProlongationLoop:;
                    }

endCycleLoop:
                    if(info.min_level == info.max_level)
                        break;

                    if(!info.is_smoother)
                    {
                        defect = Norm<vnt_l_two, true, Tag_>::value((info.d[info.max_level]));

                        *cappa = std::pow((double)(defect / initial_defect), 1.0/((Prec_)iter));

                        if (defect <= initial_defect * info.tolerance)
                            break;
                    }
                    else
                    {
                        *cappa = Prec_(-1);
                    }
                }
            }
            /*DenseVector<Prec_> result(info.x[info.max_level].size());
              copy<Tag_>(info.x[info.max_level], result);
              return result;*/
            return info.x[info.max_level];
        }

        template <typename Prec_>
        static DenseVector<Prec_> _multigrid_kernel(SparseMatrixELL<Prec_>&  system, DenseVector<Prec_>& right_hand_side, unsigned long max_levels, Prec_ * cappa, MGInfo<Prec_, SparseMatrixELL<Prec_> > & info)
        {
            bool restriction_started(false);
            // compute initial defect
            // when start vector is zero: set d = rhs (we can save one matvec here)
            // when start vector is not zero: d = rhs - A*x
            if(info.initial_zero)
            {
                info.d[info.max_level] = info.rhs[info.max_level];
            }
            else
            {
                DenseVector<Prec_> defect(Defect<Tag_>::value(info.rhs[info.max_level], info.a[info.max_level], info.x[info.max_level]));
                info.d[info.max_level] = defect;
            }

            Prec_ defect, initial_defect;
            // compute norm of initial defect, also used for relative convergence control
            if(!info.is_smoother)
            {
                defect = Norm<vnt_l_two, true, Tag_>::value(info.d[info.max_level]);
                initial_defect = defect;
            }
            else
            {
                defect = Prec_(1e8);
                initial_defect = defect;
            }
            // check if nothing needs to be done
            if(info.convergence_check && defect <= info.tolerance)
            {
                //do nothing
            }
            else
            {
                //start cycles
                unsigned long iter, current_level(info.max_level);
                unsigned long local_cycle[max_levels];
                for(iter = 1 ; iter <= info.n_max_iter ; ++iter)
                {
#ifdef SOLVER_VERBOSE
                    std::cout << iter << "th iteration" <<std::endl;
#endif
                    // set current level to the maximal level
                    current_level = info.max_level;

                    // Initialise the array providing the local cycle information
                    //V cycle!
                    for(unsigned long i(info.min_level + 1); i <= info.max_level; ++i)
                    {
                        local_cycle[i] = 1; //Using 1 to code local V cycle
                    }
                    //cycle loop
                    while(true)
                    {
                        //Restriction loop:
                        // set a flag that the restriction has just started
                        restriction_started = true;
                        while(true)
                        {
                            // for the case that we have actually only one level, do no presmoothing
                            // and jump to the coarse grid solver directly.
                            if (info.min_level == info.max_level)
                                goto endRestrictionLoop;

                            // -------------
                            // presmoothing
                            // -------------

                            if (restriction_started)
                            {
                                // When the restriction loop just started

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], (info.c[current_level]), Prec_(0.7), info.diags_inverted[current_level]);

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], (info.c[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                                Sum<Tag_>::value(info.x[current_level], info.c[current_level]);
                            }
                            else
                            {
                                // otherwise the solution process can be started directly with
                                // the cleared solution vector (as the solution vector itself represents
                                // the defect correction here)
                                //info.x[current_level] = (Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], Prec_(0.7), info.diags_inverted[current_level]);
                                //info.x[current_level] = (Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]),info.x[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                                //END NEWTEST
                            }
                            DenseVector<Prec_> defect_2(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_2;

#ifdef SOLVER_VERBOSE
                            std::cout << "-----------------------------------------------------" << std::endl;
                            std::cout << "Presmoothing ||D|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            //----------------------------------
                            //restriction ("go down one level")
                            //----------------------------------

                            --current_level;

                            // restrict defect from level icurrentLevel+1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the restricted defect vector
                            // depending on Dirichlet mask (see routine for details), and store a copy in RHS

                            /*info.d[current_level] =*/ Restriction<Tag_>::value((info.d[current_level]), (info.d[current_level + 1]), *info.macro_border_mask);
#ifdef SOLVER_VERBOSE
                            std::cout << "Restricted." << std::endl;
#endif
                            info.rhs[current_level] =(info.d[current_level]);

#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            // if we reached the coarsest level exit the restricition loop
                            if (current_level == info.min_level)
                                goto endRestrictionLoop;

                            restriction_started = false;
                        }
endRestrictionLoop:
                        // ----------------------
                        // coarse grid correction
                        // ----------------------

                        if (info.min_level == info.max_level)
                        {
                            // For the case we actually have only one MG level, only
                            // the following coarse grid correction (and no smoothing) is done.

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());

                            DenseVector<Prec_> defect_3(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_3;

                            //
                            // the MG cycle can be exited immediately
                            //
                            goto endCycleLoop;
                        }
                        else
                        {
                            // Otherwise this is a "real" coarse grid correction, which is
                            // started with a zero start vector

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());
#ifdef SOLVER_VERBOSE
                            std::cout << "Coarse Grid solver." << std::endl;
#endif
                        }

                        //-------------
                        //prolongation
                        //-------------
                        while (true)
                        {

                            //------------------------------
                            //prolongation ("go up one level")
                            //------------------------------

                            ++current_level;

                            //
                            // prolongate solution vector from level icurrentLevel-1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the prolongated correction vector
                            // depending on Dirichlet mask passed in from FEAST (see code for details)
                            //
                            /*if (current_level == 4)
                              {
                              info.c[current_level].lock(lm_write_only);
                              info.x[current_level - 1].lock(lm_read_only);
                              info.c[current_level] = Prolongation<tags::CPU>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask);
                              info.c[current_level].unlock(lm_write_only);
                              info.x[current_level - 1].unlock(lm_read_only);
                              }
                              else*/
                            {
                                /*info.c[current_level] =*/ Prolongation<Tag_, ProlType_>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask, info.prolmats[current_level]);
                            }
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongated." << std::endl;
#endif
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||c|| " << Norm<vnt_l_two, true, Tag_>::value(info.c[current_level]) << std::endl;
#endif
                            //
                            // perform adaptive coarse grid correction if required
                            //
                            Prec_ alpha;
                            /*if (info.adapt_correction_factor == 0.0)
                              {
                            //Compute dalpha = (d,c)/(Ac,c) where d is the current residual
                            //(data->d[icurrentLevel]) and c the prolongated correction (c[icurrentLevel])
                            Prec_ d1,d2;
                            d1 = DotProduct<Tag_>::value((info.d[current_level]), (info.c[current_level]));
                            DenseVector<Prec_> t(Product<Tag_>::value((info.a[current_level]), (info.c[current_level])));

                            d2 = DotProduct<Tag_>::value(t, (info.c[current_level]));

                            alpha = Prec_(d1 / d2);
                            }
                            else*/
                            {
                                alpha = info.adapt_correction_factor;
                            }

                            //
                            // add the prolongated correction to the iteration vector
                            // on the current level
                            //
                            ScaledSum<Tag_>::value((info.x[current_level]), (info.c[current_level]), alpha);

#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||x|| " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //--------------
                            //postsmoothing
                            //--------------
                            //
                            // smooth A*x = rhs based on the RHS for that level we stored during restriction

                            //(info.x[current_level]) =(Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]));
                            Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]),(info.x[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]);
                            //end NEWTEST
#ifdef SOLVER_VERBOSE
                            std::cout << "Postsmoothing ||X|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //
                            // update defect
                            //
                            DenseVector<Prec_> defect_4(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_4;
#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
                            std::cout << "-----------------------------------------------------" << std::endl;
#endif
                            // if the maximal level is reached then the MG cycle is finished,
                            // so exit the MG cycle loop
                            if (current_level == info.max_level)
                                goto endCycleLoop;

                            // If the prolongation did not reach the finest level yet, then
                            // determine how to proceed now. Dependent on the type of MG cycle the next
                            // step is a prolongation ("going up") OR a restriction ("going down").
                            if (local_cycle[current_level] == 1)
                            {
                                // local V-cycle --> go on with prolongation ("going up")
                                // this is the case when we are in a global V-cycle or at the end of
                                // a local W-cycle
                                if (true)
                                {
                                    // When we are in a global V- or F-cycle then do a V-cycle on this level
                                    // next time (in case of a global F-cycle (which does a local W-cycle first
                                    // and then a local V-cycle) it means that the local W-cycle
                                    // on this level is finished, such that next time a local V-cycle has to
                                    // be performed)
                                    local_cycle[current_level] = 1;
                                }
                                else
                                {
                                    // In case of a global W-cycle also set the local cycle to a W-cycle
                                    local_cycle[current_level] = 2;
                                }
                            }
                            else
                            {
                                // local W-cycle --> stop prolongation and do restriction again ("going down")
                                // local W-cycle becomes local V-cycle to indicate that half the W-cycle
                                // has been performed
                                local_cycle[current_level] = 1;
                                // exit prolongation loop
                                goto endProlongationLoop;
                            }
                        }
endProlongationLoop:;
                    }

endCycleLoop:
                    if(info.min_level == info.max_level)
                        break;

                    if(!info.is_smoother)
                    {
                        defect = Norm<vnt_l_two, true, Tag_>::value((info.d[info.max_level]));

                        *cappa = std::pow((double)(defect / initial_defect), 1.0/((Prec_)iter));

                        if (defect <= initial_defect * info.tolerance)
                            break;
                    }
                    else
                    {
                        *cappa = Prec_(-1);
                    }
                }
            }
            /*DenseVector<Prec_> result(info.x[info.max_level].size());
              copy<Tag_>(info.x[info.max_level], result);
              return result;*/
            return info.x[info.max_level];
        }
        public:
        template<typename Prec_>
        static inline void value(BandedMatrixQ1<Prec_>&  system,
                DenseVector<Prec_>& right_hand_side,
                DenseVector<Prec_>& x,
                unsigned long max_levels,
                Prec_ conv_rad,
                MGInfo<Prec_, BandedMatrixQ1<Prec_> > & info)
        {
            CONTEXT("When solving banded q1 linear system with MULTIGRID: ");
#ifdef SOLVER_VERBOSE_L2
            std::cout << "Calling MG solver, presmoothing=JACOBI, postsmoothing=JACOBI, coarse-grid solver=CG, datalayout=Q1, precision=FIXED" << std::endl;
#endif
            Prec_ cappa;

            //DenseVector<Prec_> initial_guess(right_hand_side.size(), Prec_(0.)); //x_0
            DenseVector<Prec_> initial_guess(x); //x_0
            DenseVector<Prec_> outer_defect(right_hand_side.size(), Prec_(0));

            // apply Dirichlet BCs for boundary nodes (semi-implicit approach)
            // note that we cleared the solution vector previously
            unsigned long sqrt_N((unsigned long)sqrt((double)right_hand_side.size()));
            for (unsigned long i(0) ; i < sqrt_N; ++i)
            {
                initial_guess[i] = right_hand_side[i];
            }
            for (unsigned long i(right_hand_side.size() - sqrt_N); i < right_hand_side.size(); ++i)
            {
                initial_guess[i] = right_hand_side[i];
            }
            for (unsigned long i(0); i < right_hand_side.size(); i += sqrt_N)
            {
                initial_guess[i] = right_hand_side[i];
            }

            /*for (unsigned long i(sqrt_N - 1); i < right_hand_side.size(); i += sqrt_N)
              {
              initial_guess[i] = right_hand_side[i];
              }*/

            info.x[info.max_level] = initial_guess;
            info.x[info.max_level] = (_multigrid_kernel<Prec_>(info.a[info.max_level], right_hand_side, max_levels, &cappa, info));
            //return info.x[info.max_level];//result;
            x = info.x[info.max_level];
        }

        template<typename Prec_>
        static inline void value(SparseMatrixELL<Prec_>&  system,
                DenseVector<Prec_>& right_hand_side,
                DenseVector<Prec_>& x,
                unsigned long max_levels,
                Prec_ conv_rad,
                MGInfo<Prec_, SparseMatrixELL<Prec_> > & info)
        {
            CONTEXT("When solving sparse ELL linear system with MULTIGRID: ");
#ifdef SOLVER_VERBOSE_L2
            std::cout << "Calling MG solver, presmoothing=JACOBI, postsmoothing=JACOBI, coarse-grid solver=CG, datalayout=ELLPACK, precision=FIXED" << std::endl;
#endif
            Prec_ cappa;

            //DenseVector<Prec_> initial_guess(right_hand_side.size(), Prec_(0.)); //x_0
            DenseVector<Prec_> initial_guess(x); //x_0
            DenseVector<Prec_> outer_defect(right_hand_side.size(), Prec_(0));

            // apply Dirichlet BCs for boundary nodes (semi-implicit approach)
            // note that we cleared the solution vector previously
            unsigned long sqrt_N((unsigned long)sqrt((double)right_hand_side.size()));
            for (unsigned long i(0) ; i < sqrt_N; ++i)
            {
                initial_guess[i] = right_hand_side[i];
            }
            for (unsigned long i(right_hand_side.size() - sqrt_N); i < right_hand_side.size(); ++i)
            {
                initial_guess[i] = right_hand_side[i];
            }
            for (unsigned long i(0); i < right_hand_side.size(); i += sqrt_N)
            {
                initial_guess[i] = right_hand_side[i];
            }

            /*for (unsigned long i(sqrt_N - 1); i < right_hand_side.size(); i += sqrt_N)
              {
              initial_guess[i] = right_hand_side[i];
              }*/

            info.x[info.max_level] = initial_guess;
            info.x[info.max_level] = (_multigrid_kernel<Prec_>(info.a[info.max_level], right_hand_side, max_levels, &cappa, info));
            //return info.x[info.max_level];//result;
            x = info.x[info.max_level];
        }
    };
    //------------------MIXED PRECISION-------------------------

    template<typename Tag_ , typename OuterTag_, typename ProlType_>
    class Multigrid<Tag_ , OuterTag_ , ProlType_, JAC, CYCLE::V, MIXED>
    {
        template <typename Prec_>
        static DenseVector<Prec_> _multigrid_kernel(unsigned long max_levels, Prec_ * cappa, MGInfo<Prec_, BandedMatrixQ1<Prec_> > & info)
        {
            bool restriction_started(false);
            // compute initial defect
            // when start vector is zero: set d = rhs (we can save one matvec here)
            // when start vector is not zero: d = rhs - A*x
            if(info.initial_zero)
            {
                info.d[info.max_level] = info.rhs[info.max_level];
            }
            else
            {
                DenseVector<Prec_> defect(Defect<Tag_>::value(info.rhs[info.max_level], info.a[info.max_level], info.x[info.max_level]));
                info.d[info.max_level] = defect;

            }

            Prec_ defect, initial_defect;
            // compute norm of initial defect, also used for relative convergence control
            if(!info.is_smoother)
            {
                defect = Norm<vnt_l_two, true, Tag_>::value(info.d[info.max_level]);
                initial_defect = defect;
            }
            else
            {
                defect = Prec_(1e8);
                initial_defect = defect;
            }
#ifdef SOLVER_VERBOSE
            std::cout << defect << std::endl;
#endif
            // check if nothing needs to be done
            if(info.convergence_check && defect <= info.tolerance)
            {
                //do nothing
            }
            else
            {
                //start cycles
                unsigned long iter, current_level(info.max_level);
                unsigned long local_cycle[max_levels];
                for(iter = 1 ; iter <= info.n_max_iter ; ++iter)
                {
#ifdef SOLVER_VERBOSE
                    std::cout << iter << "th iteration" <<std::endl;
#endif
                    // set current level to the maximal level
                    current_level = info.max_level;

                    // Initialise the array providing the local cycle information
                    //V cycle!
                    for(unsigned long i(info.min_level + 1); i <= info.max_level; ++i)
                    {
                        local_cycle[i] = 1; //Using 1 to code local V cycle
                    }
                    //cycle loop
                    while(true)
                    {
                        //Restriction loop:
                        // set a flag that the restriction has just started
                        restriction_started = true;
                        while(true)
                        {
                            // for the case that we have actually only one level, do no presmoothing
                            // and jump to the coarse grid solver directly.
                            if (info.min_level == info.max_level)
                                goto endRestrictionLoop;

                            // -------------
                            // presmoothing
                            // -------------

                            if (restriction_started)
                            {
                                // When the restriction loop just started

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], (info.c[current_level]), Prec_(0.7), info.diags_inverted[current_level]);
                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], (info.c[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);

                                Sum<Tag_>::value(info.x[current_level], info.c[current_level]);
                            }
                            else
                            {
                                // otherwise the solution process can be started directly with
                                // the cleared solution vector (as the solution vector itself represents
                                // the defect correction here)
                                //info.x[current_level] = (Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], Prec_(0.7), info.diags_inverted[current_level]);
                                //info.x[current_level] = (Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.x[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                            }
                            DenseVector<Prec_> defect_2(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_2;

#ifdef SOLVER_VERBOSE
                            std::cout << "-----------------------------------------------------" << std::endl;
                            std::cout << "Presmoothing ||D|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            //----------------------------------
                            //restriction ("go down one level")
                            //----------------------------------

                            --current_level;

                            // restrict defect from level icurrentLevel+1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the restricted defect vector
                            // depending on Dirichlet mask (see routine for details), and store a copy in RHS

                            info.d[current_level] = Restriction<Tag_>::value((info.d[current_level]), (info.d[current_level + 1]), *info.macro_border_mask);
#ifdef SOLVER_VERBOSE
                            std::cout << "Restricted." << std::endl;
#endif
                            info.rhs[current_level] =(info.d[current_level]);

#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            // if we reached the coarsest level exit the restricition loop
                            if (current_level == info.min_level)
                                goto endRestrictionLoop;

                            restriction_started = false;
                        }
endRestrictionLoop:
                        // ----------------------
                        // coarse grid correction
                        // ----------------------

                        if (info.min_level == info.max_level)
                        {
                            // For the case we actually have only one MG level, only
                            // the following coarse grid correction (and no smoothing) is done.

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());

                            DenseVector<Prec_> defect_3(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_3;

                            //
                            // the MG cycle can be exited immediately
                            //
                            goto endCycleLoop;
                        }
                        else
                        {
                            // Otherwise this is a "real" coarse grid correction, which is
                            // started with a zero start vector

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());
#ifdef SOLVER_VERBOSE
                            std::cout << "Coarse Grid solver." << std::endl;
#endif
                        }

                        //-------------
                        //prolongation
                        //-------------
                        while (true)
                        {

                            //------------------------------
                            //prolongation ("go up one level")
                            //------------------------------

                            ++current_level;

                            //
                            // prolongate solution vector from level icurrentLevel-1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the prolongated correction vector
                            // depending on Dirichlet mask passed in from FEAST (see code for details)
                            //
                            /*if (current_level == 4)
                              {
                              info.c[current_level].lock(lm_write_only);
                              info.x[current_level - 1].lock(lm_read_only);
                              info.c[current_level] = Prolongation<tags::CPU>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask);
                              info.c[current_level].unlock(lm_write_only);
                              info.x[current_level - 1].unlock(lm_read_only);
                              }
                              else*/
                            {
                                info.c[current_level] = Prolongation<Tag_, ProlType_>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask, info.prolmats[current_level]);
                            }
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongated." << std::endl;
#endif
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||c|| " << Norm<vnt_l_two, true, Tag_>::value(info.c[current_level]) << std::endl;
#endif
                            //
                            // perform adaptive coarse grid correction if required
                            //
                            Prec_ alpha;
                            if (info.adapt_correction_factor == 0.0)
                            {
                                //Compute dalpha = (d,c)/(Ac,c) where d is the current residual
                                //(data->d[icurrentLevel]) and c the prolongated correction (c[icurrentLevel])
                                Prec_ d1,d2;
                                d1 = DotProduct<Tag_>::value((info.d[current_level]), (info.c[current_level]));
                                DenseVector<Prec_> t(Product<Tag_>::value((info.a[current_level]), (info.c[current_level])));

                                d2 = DotProduct<Tag_>::value(t, (info.c[current_level]));

                                alpha = Prec_(d1 / d2);
                            }
                            else
                            {
                                alpha = info.adapt_correction_factor;
                            }

                            //
                            // add the prolongated correction to the iteration vector
                            // on the current level
                            //
                            ScaledSum<Tag_>::value((info.x[current_level]), (info.c[current_level]), alpha);

#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||x|| " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //--------------
                            //postsmoothing
                            //--------------
                            //
                            // smooth A*x = rhs based on the RHS for that level we stored during restriction
                            //
                            //(info.x[current_level]) =(Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]));
                            Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), (info.x[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]);
#ifdef SOLVER_VERBOSE
                            std::cout << "Postsmoothing ||X|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //
                            // update defect
                            //
                            DenseVector<Prec_> defect_4(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_4;
#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
                            std::cout << "-----------------------------------------------------" << std::endl;
#endif
                            // if the maximal level is reached then the MG cycle is finished,
                            // so exit the MG cycle loop
                            if (current_level == info.max_level)
                                goto endCycleLoop;

                            // If the prolongation did not reach the finest level yet, then
                            // determine how to proceed now. Dependent on the type of MG cycle the next
                            // step is a prolongation ("going up") OR a restriction ("going down").
                            if (local_cycle[current_level] == 1)
                            {
                                // local V-cycle --> go on with prolongation ("going up")
                                // this is the case when we are in a global V-cycle or at the end of
                                // a local W-cycle
                                if (true)
                                {
                                    // When we are in a global V- or F-cycle then do a V-cycle on this level
                                    // next time (in case of a global F-cycle (which does a local W-cycle first
                                    // and then a local V-cycle) it means that the local W-cycle
                                    // on this level is finished, such that next time a local V-cycle has to
                                    // be performed)
                                    local_cycle[current_level] = 1;
                                }
                                else
                                {
                                    // In case of a global W-cycle also set the local cycle to a W-cycle
                                    local_cycle[current_level] = 2;
                                }
                            }
                            else
                            {
                                // local W-cycle --> stop prolongation and do restriction again ("going down")
                                // local W-cycle becomes local V-cycle to indicate that half the W-cycle
                                // has been performed
                                local_cycle[current_level] = 1;
                                // exit prolongation loop
                                goto endProlongationLoop;
                            }
                        }
endProlongationLoop:;
                    }

endCycleLoop:
                    if(info.min_level == info.max_level)
                        break;

                    if(!info.is_smoother)
                    {
                        defect = Norm<vnt_l_two, true, Tag_>::value((info.d[info.max_level]));

                        *cappa = std::pow((double)(defect / initial_defect), 1.0/((Prec_)iter));

                        if (defect <= initial_defect * info.tolerance)
                            break;
                    }
                    else
                    {
                        *cappa = Prec_(-1);
                    }
                }
                // write number of iterations needed until convergence
                if (iter > info.n_max_iter)
                    info.inner_iterations = iter - 1;
                else
                    info.inner_iterations = iter;
            }
            DenseVector<Prec_> result(info.x[info.max_level].size());
            copy<OuterTag_>(info.x[info.max_level], result);

            return result;
        }

        template <typename Prec_>
        static DenseVector<Prec_> _multigrid_kernel(unsigned long max_levels, Prec_ * cappa, MGInfo<Prec_, SparseMatrixELL<Prec_> > & info)
        {
            bool restriction_started(false);
            // compute initial defect
            // when start vector is zero: set d = rhs (we can save one matvec here)
            // when start vector is not zero: d = rhs - A*x
            if(info.initial_zero)
            {
                info.d[info.max_level] = info.rhs[info.max_level];
            }
            else
            {
                DenseVector<Prec_> defect(Defect<Tag_>::value(info.rhs[info.max_level], info.a[info.max_level], info.x[info.max_level]));
                info.d[info.max_level] = defect;

            }

            Prec_ defect, initial_defect;
            // compute norm of initial defect, also used for relative convergence control
            if(!info.is_smoother)
            {
                defect = Norm<vnt_l_two, true, Tag_>::value(info.d[info.max_level]);
                initial_defect = defect;
            }
            else
            {
                defect = Prec_(1e8);
                initial_defect = defect;
            }
#ifdef SOLVER_VERBOSE
            std::cout << defect << std::endl;
#endif
            // check if nothing needs to be done
            if(info.convergence_check && defect <= info.tolerance)
            {
                //do nothing
            }
            else
            {
                //start cycles
                unsigned long iter, current_level(info.max_level);
                unsigned long local_cycle[max_levels];
                for(iter = 1 ; iter <= info.n_max_iter ; ++iter)
                {
#ifdef SOLVER_VERBOSE
                    std::cout << iter << "th iteration" <<std::endl;
#endif
                    // set current level to the maximal level
                    current_level = info.max_level;

                    // Initialise the array providing the local cycle information
                    //V cycle!
                    for(unsigned long i(info.min_level + 1); i <= info.max_level; ++i)
                    {
                        local_cycle[i] = 1; //Using 1 to code local V cycle
                    }
                    //cycle loop
                    while(true)
                    {
                        //Restriction loop:
                        // set a flag that the restriction has just started
                        restriction_started = true;
                        while(true)
                        {
                            // for the case that we have actually only one level, do no presmoothing
                            // and jump to the coarse grid solver directly.
                            if (info.min_level == info.max_level)
                                goto endRestrictionLoop;

                            // -------------
                            // presmoothing
                            // -------------

                            if (restriction_started)
                            {
                                // When the restriction loop just started

                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], (info.c[current_level]), Prec_(0.7), info.diags_inverted[current_level]);
                                //(info.c[current_level]) = (Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.c[current_level] , info.a[current_level], info.d[current_level], (info.c[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);

                                Sum<Tag_>::value(info.x[current_level], info.c[current_level]);
                            }
                            else
                            {
                                // otherwise the solution process can be started directly with
                                // the cleared solution vector (as the solution vector itself represents
                                // the defect correction here)
                                //info.x[current_level] = (Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], Prec_(0.7), info.diags_inverted[current_level]);
                                //info.x[current_level] = (Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]));
                                Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.d[current_level]), info.x[current_level], info.n_pre_smooth - 1, Prec_(0.7), info.diags_inverted[current_level]);
                            }
                            DenseVector<Prec_> defect_2(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_2;

#ifdef SOLVER_VERBOSE
                            std::cout << "-----------------------------------------------------" << std::endl;
                            std::cout << "Presmoothing ||D|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            //----------------------------------
                            //restriction ("go down one level")
                            //----------------------------------

                            --current_level;

                            // restrict defect from level icurrentLevel+1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the restricted defect vector
                            // depending on Dirichlet mask (see routine for details), and store a copy in RHS

                            info.d[current_level] = Restriction<Tag_>::value((info.d[current_level]), (info.d[current_level + 1]), *info.macro_border_mask);
#ifdef SOLVER_VERBOSE
                            std::cout << "Restricted." << std::endl;
#endif
                            info.rhs[current_level] =(info.d[current_level]);

#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
#endif
                            // if we reached the coarsest level exit the restricition loop
                            if (current_level == info.min_level)
                                goto endRestrictionLoop;

                            restriction_started = false;
                        }
endRestrictionLoop:
                        // ----------------------
                        // coarse grid correction
                        // ----------------------

                        if (info.min_level == info.max_level)
                        {
                            // For the case we actually have only one MG level, only
                            // the following coarse grid correction (and no smoothing) is done.

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());

                            DenseVector<Prec_> defect_3(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_3;

                            //
                            // the MG cycle can be exited immediately
                            //
                            goto endCycleLoop;
                        }
                        else
                        {
                            // Otherwise this is a "real" coarse grid correction, which is
                            // started with a zero start vector

                            //(info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
                            ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), info.x[current_level], std::numeric_limits<Prec_>::epsilon());
#ifdef SOLVER_VERBOSE
                            std::cout << "Coarse Grid solver." << std::endl;
#endif
                        }

                        //-------------
                        //prolongation
                        //-------------
                        while (true)
                        {

                            //------------------------------
                            //prolongation ("go up one level")
                            //------------------------------

                            ++current_level;

                            //
                            // prolongate solution vector from level icurrentLevel-1 to icurrentLevel,
                            // set homogeneous Dirichlet boundary conditions in the prolongated correction vector
                            // depending on Dirichlet mask passed in from FEAST (see code for details)
                            //
                            /*if (current_level == 4)
                              {
                              info.c[current_level].lock(lm_write_only);
                              info.x[current_level - 1].lock(lm_read_only);
                              info.c[current_level] = Prolongation<tags::CPU>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask);
                              info.c[current_level].unlock(lm_write_only);
                              info.x[current_level - 1].unlock(lm_read_only);
                              }
                              else*/
                            {
                                info.c[current_level] = Prolongation<Tag_, ProlType_>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask, info.prolmats[current_level]);
                            }
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongated." << std::endl;
#endif
#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||c|| " << Norm<vnt_l_two, true, Tag_>::value(info.c[current_level]) << std::endl;
#endif
                            //
                            // perform adaptive coarse grid correction if required
                            //
                            Prec_ alpha;
                            if (info.adapt_correction_factor == 0.0)
                            {
                                //Compute dalpha = (d,c)/(Ac,c) where d is the current residual
                                //(data->d[icurrentLevel]) and c the prolongated correction (c[icurrentLevel])
                                Prec_ d1,d2;
                                d1 = DotProduct<Tag_>::value((info.d[current_level]), (info.c[current_level]));
                                DenseVector<Prec_> t(info.c[current_level].size(), Prec_(0));
                                Product<Tag_>::value(t , (info.a[current_level]), (info.c[current_level]));

                                d2 = DotProduct<Tag_>::value(t, (info.c[current_level]));

                                alpha = Prec_(d1 / d2);
                            }
                            else
                            {
                                alpha = info.adapt_correction_factor;
                            }

                            //
                            // add the prolongated correction to the iteration vector
                            // on the current level
                            //
                            ScaledSum<Tag_>::value((info.x[current_level]), (info.c[current_level]), alpha);

#ifdef SOLVER_VERBOSE
                            std::cout << "Prolongation on level " << current_level << " ||x|| " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //--------------
                            //postsmoothing
                            //--------------
                            //
                            // smooth A*x = rhs based on the RHS for that level we stored during restriction
                            //
                            //(info.x[current_level]) =(Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]));
                            Jacobi<Tag_>::value(info.x[current_level], (info.a[current_level]), (info.rhs[current_level]), (info.x[current_level]), info.n_pre_smooth, Prec_(0.7), info.diags_inverted[current_level]);
#ifdef SOLVER_VERBOSE
                            std::cout << "Postsmoothing ||X|| on level " << current_level << " " << Norm<vnt_l_two, true, Tag_>::value(info.x[current_level]) << std::endl;
#endif
                            //
                            // update defect
                            //
                            DenseVector<Prec_> defect_4(Defect<Tag_>::value(info.rhs[current_level], info.a[current_level], info.x[current_level]));
                            info.d[current_level] = defect_4;
#ifdef SOLVER_VERBOSE
                            std::cout << "Defect on level " << current_level << "||D||: " << Norm<vnt_l_two, true, Tag_>::value(info.d[current_level]) << std::endl;
                            std::cout << "-----------------------------------------------------" << std::endl;
#endif
                            // if the maximal level is reached then the MG cycle is finished,
                            // so exit the MG cycle loop
                            if (current_level == info.max_level)
                                goto endCycleLoop;

                            // If the prolongation did not reach the finest level yet, then
                            // determine how to proceed now. Dependent on the type of MG cycle the next
                            // step is a prolongation ("going up") OR a restriction ("going down").
                            if (local_cycle[current_level] == 1)
                            {
                                // local V-cycle --> go on with prolongation ("going up")
                                // this is the case when we are in a global V-cycle or at the end of
                                // a local W-cycle
                                if (true)
                                {
                                    // When we are in a global V- or F-cycle then do a V-cycle on this level
                                    // next time (in case of a global F-cycle (which does a local W-cycle first
                                    // and then a local V-cycle) it means that the local W-cycle
                                    // on this level is finished, such that next time a local V-cycle has to
                                    // be performed)
                                    local_cycle[current_level] = 1;
                                }
                                else
                                {
                                    // In case of a global W-cycle also set the local cycle to a W-cycle
                                    local_cycle[current_level] = 2;
                                }
                            }
                            else
                            {
                                // local W-cycle --> stop prolongation and do restriction again ("going down")
                                // local W-cycle becomes local V-cycle to indicate that half the W-cycle
                                // has been performed
                                local_cycle[current_level] = 1;
                                // exit prolongation loop
                                goto endProlongationLoop;
                            }
                        }
endProlongationLoop:;
                    }

endCycleLoop:
                    if(info.min_level == info.max_level)
                        break;

                    if(!info.is_smoother)
                    {
                        defect = Norm<vnt_l_two, true, Tag_>::value((info.d[info.max_level]));

                        *cappa = std::pow((double)(defect / initial_defect), 1.0/((Prec_)iter));

                        if (defect <= initial_defect * info.tolerance)
                            break;
                    }
                    else
                    {
                        *cappa = Prec_(-1);
                    }
                }
                // write number of iterations needed until convergence
                if (iter > info.n_max_iter)
                    info.inner_iterations = iter - 1;
                else
                    info.inner_iterations = iter;
            }
            DenseVector<Prec_> result(info.x[info.max_level].size());
            copy<OuterTag_>(info.x[info.max_level], result);

            return result;
        }

        public:
        template<typename InnerPrec_, typename OuterPrec_>
        static inline void value(BandedMatrixQ1<OuterPrec_>&  system,
                DenseVector<OuterPrec_>& right_hand_side,
                DenseVector<OuterPrec_>& x,
                unsigned long max_levels,
                OuterPrec_ conv_rad,
                MGInfo<InnerPrec_, BandedMatrixQ1<InnerPrec_> > & info)
            {
                CONTEXT("When solving banded q1 linear system with MULTIGRID: ");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling MG solver, presmoothing=JACOBI, postsmoothing=JACOBI, coarse-grid solver=CG, datalayout=Q1, precision=MIXED" << std::endl;
#endif
#ifdef SOLVER_BENCHMARK
                TimeStamp ab, ae;
                ab.take();
                TimeStamp pb, pe;
                pb.take();
#endif
                InnerPrec_ cappa;

                // current and initial defect
                OuterPrec_ initial_defect;
                DenseVector<OuterPrec_> result(right_hand_side.size(), OuterPrec_(0)); //final result
                //DenseVector<OuterPrec_> initial_guess(right_hand_side.size(), OuterPrec_(0)); //x_0
                DenseVector<OuterPrec_> initial_guess(x); //x_0

                // apply Dirichlet BCs for boundary nodes (semi-implicit approach)
                // note that we cleared the solution vector previously
                unsigned long sqrt_N((unsigned long)sqrt((double)right_hand_side.size()));
                for (unsigned long i(0) ; i < sqrt_N; ++i)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                for (unsigned long i(right_hand_side.size() - sqrt_N); i < right_hand_side.size(); ++i)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                for (unsigned long i(0); i < right_hand_side.size(); i += sqrt_N)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                /*for (unsigned long i(sqrt_N - 1); i < right_hand_side.size(); i += sqrt_N)
                  {
                  initial_guess[i] = right_hand_side[i];
                  }*/

                convert<OuterTag_>(info.x[info.max_level], initial_guess);

                unsigned long inner_iterations(0);
                unsigned long outer_iterations(0);
                OuterPrec_ scale_factor(1.0);

                // calc initial defect
                // D = B - Ax, d0 = ||D||
                DenseVector<OuterPrec_> product(Product<OuterTag_>::value(system, initial_guess));
                DenseVector<OuterPrec_> outer_defect(right_hand_side.size());
                copy<OuterTag_>(right_hand_side, outer_defect);
                Difference<OuterTag_>::value(outer_defect, product);
                initial_defect = Norm<vnt_l_two, true, OuterTag_>::value(outer_defect);

                OuterPrec_ def_norm;
                OuterPrec_ inv;
#ifdef SOLVER_BENCHMARK
                pe.take();
                std::cout << "PreProc TOE: "<< (pe.sec() - pb.sec()) + (pe.usec() - pb.usec())/1e6 << std::endl;
#endif
                DenseVector<OuterPrec_> x_outer(result.size());
                while(inner_iterations < 8)
                {
#ifdef SOLVER_BENCHMARK
                    TimeStamp ob, oe, ib, ie;
                    ob.take();
#endif
                    // set defect as RHS to inner solver
                    convert<OuterTag_>(info.rhs[info.max_level], outer_defect);

                    // run inner solver as long as neccessary

#ifdef SOLVER_VERBOSE
                    std::cout << inner_iterations << "th iteration (outer)!" << std::endl;
#endif
#ifdef SOLVER_BENCHMARK
                    ib.take();
#endif
                    info.x[info.max_level] = (_multigrid_kernel<InnerPrec_>(max_levels, &cappa, info));
#ifdef SOLVER_BENCHMARK
                    ie.take();
#endif

                    inner_iterations += info.inner_iterations;

                    // get "solution" and update outer solution

                    convert<OuterTag_>(x_outer, info.x[info.max_level]);
                    ScaledSum<OuterTag_>::value(result, x_outer, scale_factor);


                    // calculate defect
                    outer_defect = Defect<OuterTag_>::value(right_hand_side, system, result);

                    // calc norm of defect
                    def_norm = Norm<vnt_l_two, true, OuterTag_>::value(outer_defect);
#ifdef SOLVER_VERBOSE
                    std::cout << "defnorm: " << def_norm << std::endl;
                    std::cout << "initial_defect: " << initial_defect << std::endl;
#endif
                    // check for convergence
                    OuterPrec_ outer_eps(1e-8);
                    if (def_norm  < outer_eps * initial_defect)
                    {
#ifdef SOLVER_BENCHMARK
                        oe.take();
                        std::cout << "Outer TOE: "<< (oe.sec() - ob.sec()) + (oe.usec() - ob.usec())/1e6<< std::endl;
                        std::cout << "Inner TOE: "<< (ie.sec() - ib.sec()) + (ie.usec() - ib.usec())/1e6<< std::endl;
#endif
                        break;
                    }
                    // scale defect
                    scale_factor = def_norm;

                    inv = (1.0 / scale_factor);
                    Scale<OuterTag_>::value(outer_defect, inv);

                    outer_iterations++;
#ifdef SOLVER_BENCHMARK
                    oe.take();
                    std::cout << "Outer TOE: "<< (oe.sec() - ob.sec()) + (oe.usec() - ob.usec())/1e6<< std::endl;
                    std::cout << "Inner TOE: "<< (ie.sec() - ib.sec()) + (ie.usec() - ib.usec())/1e6<< std::endl;
#endif
                }
                //#ifdef SOLVER_VERBOSE
                std::cout << "TN of outer iters: " << outer_iterations << std::endl;
                std::cout << "TN of inner iters: " << inner_iterations << std::endl;
                //#endif
#ifdef SOLVER_BENCHMARK
                ae.take();
                std::cout << "All TOE: "<< (ae.sec() - ab.sec()) + (ae.usec() - ab.usec())/1e6 << std::endl;
#endif
                //return result;
                x = result;
            }

            template<typename InnerPrec_, typename OuterPrec_>
            static inline void value(SparseMatrixELL<OuterPrec_>&  system,
                DenseVector<OuterPrec_>& right_hand_side,
                DenseVector<OuterPrec_>& x,
                unsigned long max_levels,
                OuterPrec_ conv_rad,
                MGInfo<InnerPrec_, SparseMatrixELL<InnerPrec_> > & info)
            {
                CONTEXT("When solving banded q1 linear system with MULTIGRID: ");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling MG solver, presmoothing=JACOBI, postsmoothing=JACOBI, coarse-grid solver=CG, datalayout=ELL, precision=MIXED" << std::endl;
#endif
#ifdef SOLVER_BENCHMARK
                TimeStamp ab, ae;
                ab.take();
                TimeStamp pb, pe;
                pb.take();
#endif
                InnerPrec_ cappa;

                // current and initial defect
                OuterPrec_ initial_defect;
                DenseVector<OuterPrec_> result(right_hand_side.size(), OuterPrec_(0)); //final result
                //DenseVector<OuterPrec_> initial_guess(right_hand_side.size(), OuterPrec_(0)); //x_0
                DenseVector<OuterPrec_> initial_guess(x); //x_0

                // apply Dirichlet BCs for boundary nodes (semi-implicit approach)
                // note that we cleared the solution vector previously
                unsigned long sqrt_N((unsigned long)sqrt((double)right_hand_side.size()));
                for (unsigned long i(0) ; i < sqrt_N; ++i)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                for (unsigned long i(right_hand_side.size() - sqrt_N); i < right_hand_side.size(); ++i)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                for (unsigned long i(0); i < right_hand_side.size(); i += sqrt_N)
                {
                    initial_guess[i] = right_hand_side[i];
                }
                /*for (unsigned long i(sqrt_N - 1); i < right_hand_side.size(); i += sqrt_N)
                  {
                  initial_guess[i] = right_hand_side[i];
                  }*/

                convert<OuterTag_>(info.x[info.max_level], initial_guess);

                unsigned long inner_iterations(0);
                unsigned long outer_iterations(0);
                OuterPrec_ scale_factor(1.0);

                // calc initial defect
                // D = B - Ax, d0 = ||D||
                DenseVector<OuterPrec_> product(right_hand_side.size(), OuterPrec_(0));
                Product<OuterTag_>::value(product, system, initial_guess);
                DenseVector<OuterPrec_> outer_defect(right_hand_side.size());
                copy<OuterTag_>(right_hand_side, outer_defect);
                Difference<OuterTag_>::value(outer_defect, product);
                initial_defect = Norm<vnt_l_two, true, OuterTag_>::value(outer_defect);

                OuterPrec_ def_norm;
                OuterPrec_ inv;
#ifdef SOLVER_BENCHMARK
                pe.take();
                std::cout << "PreProc TOE: "<< (pe.sec() - pb.sec()) + (pe.usec() - pb.usec())/1e6 << std::endl;
#endif
                DenseVector<OuterPrec_> x_outer(result.size());
                while(inner_iterations < 8)
                {
#ifdef SOLVER_BENCHMARK
                    TimeStamp ob, oe, ib, ie;
                    ob.take();
#endif
                    // set defect as RHS to inner solver
                    convert<OuterTag_>(info.rhs[info.max_level], outer_defect);

                    // run inner solver as long as neccessary

#ifdef SOLVER_VERBOSE
                    std::cout << inner_iterations << "th iteration (outer)!" << std::endl;
#endif
#ifdef SOLVER_BENCHMARK
                    ib.take();
#endif
                    info.x[info.max_level] = (_multigrid_kernel<InnerPrec_>(max_levels, &cappa, info));
#ifdef SOLVER_BENCHMARK
                    ie.take();
#endif

                    inner_iterations += info.inner_iterations;

                    // get "solution" and update outer solution

                    convert<OuterTag_>(x_outer, info.x[info.max_level]);
                    ScaledSum<OuterTag_>::value(result, x_outer, scale_factor);


                    // calculate defect
                    outer_defect = Defect<OuterTag_>::value(right_hand_side, system, result);

                    // calc norm of defect
                    def_norm = Norm<vnt_l_two, true, OuterTag_>::value(outer_defect);
#ifdef SOLVER_VERBOSE
                    std::cout << "defnorm: " << def_norm << std::endl;
                    std::cout << "initial_defect: " << initial_defect << std::endl;
#endif
                    // check for convergence
                    OuterPrec_ outer_eps(1e-8);
                    if (def_norm  < outer_eps * initial_defect)
                    {
#ifdef SOLVER_BENCHMARK
                        oe.take();
                        std::cout << "Outer TOE: "<< (oe.sec() - ob.sec()) + (oe.usec() - ob.usec())/1e6<< std::endl;
                        std::cout << "Inner TOE: "<< (ie.sec() - ib.sec()) + (ie.usec() - ib.usec())/1e6<< std::endl;
#endif
                        break;
                    }
                    // scale defect
                    scale_factor = def_norm;

                    inv = (1.0 / scale_factor);
                    Scale<OuterTag_>::value(outer_defect, inv);

                    outer_iterations++;
#ifdef SOLVER_BENCHMARK
                    oe.take();
                    std::cout << "Outer TOE: "<< (oe.sec() - ob.sec()) + (oe.usec() - ob.usec())/1e6<< std::endl;
                    std::cout << "Inner TOE: "<< (ie.sec() - ib.sec()) + (ie.usec() - ib.usec())/1e6<< std::endl;
#endif
                }
                //#ifdef SOLVER_VERBOSE
                std::cout << "TN of outer iters: " << outer_iterations << std::endl;
                std::cout << "TN of inner iters: " << inner_iterations << std::endl;
                //#endif
#ifdef SOLVER_BENCHMARK
                ae.take();
                std::cout << "All TOE: "<< (ae.sec() - ab.sec()) + (ae.usec() - ab.usec())/1e6 << std::endl;
#endif
                //return result;
                x = result;
            }
    };
}
#endif
