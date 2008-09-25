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

//#define SOLVER_VERBOSE 1
using namespace methods;
namespace honei
{

    template <typename Prec_>
    struct MGInfo
    {
        public:
            //configuration constants:
            bool is_smoother;
            DenseVector<unsigned long>* macro_border_mask;
            unsigned long min_level;
            unsigned long max_level;

            unsigned long  n_max_iter;
            bool initial_zero;
            Prec_ tolerance;
            bool convergence_check;

            unsigned long  n_pre_smooth;
            unsigned long  n_post_smooth;
            unsigned long  n_max_iter_coarse;
            Prec_ tolerance_coarse;
            Prec_ adapt_correction_factor;
            unsigned long n_truncate_d_o_f;
            unsigned long n_threshold_d_o_f;

            std::vector<DenseVector<Prec_> > c;
            std::vector<DenseVector<Prec_> > d;
            std::vector<DenseVector<Prec_> > rhs;
            std::vector<DenseVector<Prec_> > x;

            std::vector<BandedMatrixQ1<Prec_> > a;

    };

    template<typename Tag_, typename SmootherType_, typename CycleType_, typename Mode_>
        struct Multigrid
        {
        };

    template<typename Tag_>
        class Multigrid<Tag_, JAC, CYCLE::V, FIXED>
        {
                template <typename Prec_>
                    static DenseVector<Prec_> _multigrid_kernel(BandedMatrixQ1<Prec_>&  system, DenseVector<Prec_>& right_hand_side, unsigned long max_levels, Prec_ * cappa, MGInfo<Prec_> & info)
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

                                            (info.c[current_level]) = (Jacobi<Tag_>::value(info.a[current_level], info.d[current_level], Prec_(0.7)));
                                            //DenseVector<Prec_> null(info.x[current_level].size() , Prec_(0));
                                            DenseVector<Prec_> temp_jac(info.c[current_level].size());
                                            copy<Tag_>(info.c[current_level], temp_jac);
                                            (info.c[current_level]) = (Jacobi<Tag_>::value(temp_jac , info.a[current_level], info.d[current_level], info.n_pre_smooth - 1, Prec_(0.7)));

                                            Sum<Tag_>::value(info.x[current_level], info.c[current_level]);
                                        }
                                        else
                                        {
                                            // otherwise the solution process can be started directly with
                                            // the cleared solution vector (as the solution vector itself represents
                                            // the defect correction here)
                                            info.x[current_level] = (Jacobi<Tag_>::value((info.a[current_level]), (info.d[current_level]), Prec_(0.7)));
                                            //DenseVector<Prec_> null(info.x[current_level].size() , Prec_(0));
                                            DenseVector<Prec_> temp_jac(info.x[current_level].size());
                                            copy<Tag_>(info.x[current_level], temp_jac);
                                            info.x[current_level] = (Jacobi<Tag_>::value(temp_jac, (info.a[current_level]), (info.d[current_level]), info.n_pre_smooth - 1, Prec_(0.7)));
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

                                        (info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));

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

                                        (info.x[current_level]) =(ConjugateGradients<Tag_, NONE>::value((info.a[current_level]), (info.d[current_level]), std::numeric_limits<Prec_>::epsilon()));
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
                                            info.c[current_level] = Prolongation<Tag_>::value((info.c[current_level]), (info.x[current_level - 1]), *info.macro_border_mask);
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
                                        DenseVector<Prec_> temp_jac(info.x[current_level].size());
                                        copy<Tag_>(info.x[current_level], temp_jac);
                                        (info.x[current_level]) =(Jacobi<Tag_>::value(temp_jac, (info.a[current_level]), (info.rhs[current_level]), info.n_pre_smooth, Prec_(0.7)));
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

                                    *cappa = pow(defect / initial_defect, 1.0/((Prec_)iter));

                                    if (defect <= initial_defect * info.tolerance)
                                        break;
                                }
                                else
                                {
                                    *cappa = Prec_(-1);
                                }
                            }
                        }
                        DenseVector<Prec_> result(info.x[info.max_level].size());
                        copy<Tag_>(info.x[info.max_level], result);
                        return result;
                    }

            public:
                template<typename Prec_>
                    static DenseVector<Prec_> value(BandedMatrixQ1<Prec_>&  system, DenseVector<Prec_>& right_hand_side, unsigned long max_levels, Prec_ conv_rad, MGInfo<Prec_> & info)
                    {
                        Prec_ cappa;

                        // cycle control
                        unsigned long iter;
                        bool restriction_started;
                        unsigned long local_cycle[max_levels];
                        for(unsigned long i(0); i < max_levels; ++i)
                        {
                            local_cycle[i] = 1;
                        }

                        // current and initial defect
                        Prec_ defect, initial_defect;
                        DenseVector<Prec_> result(right_hand_side.size(), Prec_(0)); //final result


                        DenseVector<Prec_> initial_guess(right_hand_side.size(), Prec_(0)); //x_0
                        DenseVector<Prec_> outer_defect(right_hand_side.size(), Prec_(0));
                        DenseVector<Prec_> temp_vector(right_hand_side.size(), Prec_(0)); //delete if unneeded

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

                        for (unsigned long i(sqrt_N - 1); i < right_hand_side.size(); i += sqrt_N)
                        {
                            initial_guess[i] = right_hand_side[i];
                        }

                        unsigned long timing_loop(1);
                        for(unsigned long timing(0); timing < timing_loop; ++timing)
                        {
                            unsigned long inner_iterations(0);
                            unsigned long outer_iterations(1);
                            Prec_ scale_factor(1.0);

                            // calc initial defect
                            // D = B - Ax, d0 = ||D||
                            DenseVector<Prec_> product(Product<Tag_>::value(info.a[info.max_level], initial_guess));
                            DenseVector<Prec_> rhs_c(right_hand_side.size());
                            copy<Tag_>(right_hand_side, rhs_c);
                            Difference<Tag_>::value(rhs_c, product);
                            outer_defect = rhs_c;
                            initial_defect = Norm<vnt_l_two, true, Tag_>::value(outer_defect);

                            Prec_ def_norm;
                            Prec_ inv;
                            unsigned long step_iterations;
                            while(inner_iterations < 1)
                            {
                                // set defect as RHS to inner solver
                                temp_vector.lock(lm_write_only);
                                outer_defect.lock(lm_read_only);
                                for (unsigned long i(0); i < right_hand_side.size(); ++i)
                                {
                                    temp_vector[i] = (Prec_)outer_defect[i];
                                }
                                temp_vector.unlock(lm_write_only);
                                outer_defect.unlock(lm_read_only);
                                copy<Tag_>(temp_vector, info.rhs[info.max_level]);

                                // run inner solver as long as neccessary
#ifdef SOLVER_VERBOSE
                                std::cout << inner_iterations << "th iteration (outer)!" << std::endl;
#endif
                                info.x[info.max_level] = (_multigrid_kernel<Prec_>(info.a[info.max_level], right_hand_side, max_levels, &cappa, info));
                                inner_iterations += 1; //Markus: for now, this is ok, later: add kernel iterations

                                // get "solution" and update outer solution
                                copy<Tag_>(info.x[info.max_level], temp_vector);

                                result.lock(lm_write_only);
                                temp_vector.lock(lm_read_only);
                                for (unsigned long i(0); i < right_hand_side.size(); ++i)
                                {
                                    result[i] += scale_factor * (Prec_)(temp_vector)[i];
                                }
                                result.unlock(lm_write_only);
                                temp_vector.unlock(lm_read_only);


                                // calculate defect

                                DenseVector<Prec_> defect_outer(Defect<Tag_>::value(right_hand_side, info.a[info.max_level], result));
                                outer_defect = defect_outer;

                                // calc norm of defect
                                def_norm = Norm<vnt_l_two, true, Tag_>::value(outer_defect);
#ifdef SOLVER_VERBOSE
                                std::cout << "defnorm: " << def_norm << std::endl;
                                std::cout << "initial_defect: " << initial_defect << std::endl;
#endif
                                // check for convergence
                                Prec_ outer_eps(1e-8);
                                if (def_norm  < outer_eps * initial_defect)
                                {
                                    break;
                                }
                                // scale defect
                                scale_factor = def_norm;

                                inv = (1.0 / scale_factor);
                                Scale<Tag_>::value(outer_defect, inv);

                                outer_iterations++;
                            }
                        }

                        return info.x[info.max_level];//result;
                    }
        };
}
#endif
