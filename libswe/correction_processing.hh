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


#ifndef LIBSWE_GUARD_CORRECTION_PROCESSING_HH
#define LIBSWE_GUARD_CORRECTION_PROCESSING_HH 1

#include <libla/dense_vector.hh>
#include <libla/product.hh>
#include <libswe/directions.hh>
#include <libswe/boundary_types.hh>

/**
 * \file
 * Implementation of correction processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    template <typename BoundaryType_, typename Tag_>
    struct CorrectionProcessing
    {
    };

    /**
     * \brief Implementation of correction processing.
     *
     * \ingroup grplibswe
     *
     **/
    template <>
    struct CorrectionProcessing<boundaries::REFLECT, tags::CPU>
    {
        public:
            /**
             * \brief Implementation of correction processing.
             *
             * \param predictedu Vector containing the prediction of u.
             * \param predictedv Vector containing the prediction of v.
             * \param predictedw Vector containing the prediction of w.
             * \param u The target vector u.
             * \param v The target vector v.
             * \param w The target vector w.
             * \param d_width Grid width.
             * \param d_height Grid height.
             * \param height The result matrix.
             *
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline void value(DenseVector<WorkPrec_>& predictedu, DenseVector<WorkPrec_>& predictedv, DenseVector<WorkPrec_>& predictedw, DenseVector<WorkPrec_>& u, DenseVector<WorkPrec_>& v, DenseVector<WorkPrec_>& w, unsigned long d_width, unsigned long d_height, DenseMatrix<WorkPrec_> & height)
            {
                CONTEXT("When processing RelaxSolver correction.");

                ///correct first 2(w+4)+2 ghost cells (tripels)
                typename DenseVector<WorkPrec_>::ElementIterator iter(u.begin_elements());
                while(iter.index()<(6*(d_width+4)+6))
                {
                    ++iter;
                }
#ifdef SOLVER_VERBOSE
                cout << stringify(iter.index()) << endl;
#endif
                unsigned long count =0;
                ///Iterate through predicted u,v,w - vectors, compute weighted sum , read out h_ij, care about ghost cells.
                typename DenseMatrix<WorkPrec_>::ElementIterator h(height.begin_elements());
                while(iter.index()<3*((d_height+2)*(d_width+4)))
                {
                    if(count < d_width)
                    {
                        WorkPrec_ a = (u)[iter.index()];
                        bool warning = false;
                        if(predictedu[iter.index()] + a < 0)
                            warning = true;

                        (u)[iter.index()] = fabs(WorkPrec_(0.5)*(predictedu[iter.index()] + a));
                        *h = (u)[iter.index()];
                        ++h;
                        ++count;

                        (v)[iter.index()] = WorkPrec_(0.5)*(predictedv[iter.index()]+ (v)[iter.index()]);
                        (w)[iter.index()] = WorkPrec_(0.5)*(predictedw[iter.index()]+ (w)[iter.index()]);
                        ++iter;
                        if(warning)
                            (u)[iter.index()] = fabs(WorkPrec_(0.5)*(predictedu[iter.index()]+ (u)[iter.index()]));
                        else
                            (u)[iter.index()] = WorkPrec_(0.5)*(predictedu[iter.index()]+ (u)[iter.index()]);

                        (v)[iter.index()] = WorkPrec_(0.5)*(predictedv[iter.index()]+ (v)[iter.index()]);
                        (w)[iter.index()] = WorkPrec_(0.5)*(predictedw[iter.index()]+ (w)[iter.index()]);
                        ++iter;
                        if(warning)
                            (u)[iter.index()] = fabs(WorkPrec_(0.5)*(predictedu[iter.index()]+ (u)[iter.index()]));
                        else
                            (u)[iter.index()] = WorkPrec_(0.5)*(predictedu[iter.index()]+ (u)[iter.index()]);

                        (v)[iter.index()] = WorkPrec_(0.5)*(predictedv[iter.index()]+ (v)[iter.index()]);
                        (w)[iter.index()] = WorkPrec_(0.5)*(predictedw[iter.index()]+ (w)[iter.index()]);
                        ++iter;
                    }
                    else
                    {
                        for (unsigned long k=0; k<12; ++k)
                        {
                            ++iter;
                        }
                        count = 0;
                    }
#ifdef SOLVER_VERBOSE
                    cout << stringify(count)<<endl;
#endif
                }
#ifdef SOLVER_VERBOSE
                std::cout << "Finished Correction.\n";
#endif

            }
    };
}
#endif
