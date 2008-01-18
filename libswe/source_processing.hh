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


#ifndef LIBSWE_GUARD_SOURCE_PROCESSING_HH
#define LIBSWE_GUARD_SOURCE_PROCESSING_HH 1

#include <libla/dense_vector.hh>
#include <libla/product.hh>
#include <libswe/directions.hh>
/**
 * \file
 * Implementation of source processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    namespace source_types
    {
        class SIMPLE
        {
        };

        class INTEGRAL
        {
        };
    }

    template <typename Type_, typename Tag_>
    struct SourceProcessing
    {
    };

    /**
     * \brief Implementation of source processing.
     *
     * \ingroup grplibswe
     *
     **/
    template <typename Tag_>
    struct SourceProcessing<source_types::SIMPLE, Tag_>
    {
        private:
            /**
             * \brief Implementation of elementwise source processing.
             *
             * \param h The height at current position.
             * \param q1 The height times velocity in x direction.
             * \param q2 The height times velocity in y direction.
             * \param slope_x The slope value in x direction.
             * \param slope_y The slope value in y direction.
             * \param manning_n_squared The manning constant squared.
             *
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> _source(WorkPrec_ h, WorkPrec_ q1, WorkPrec_ q2, WorkPrec_ slope_x, WorkPrec_ slope_y, WorkPrec_ manning_n_squared)
            {
                CONTEXT("When processing RelaxSolver source, type: simple");

                DenseVector<WorkPrec_> result((unsigned long)(3), WorkPrec_(0));

                if (fabs(h) >= std::numeric_limits<WorkPrec_>::epsilon())
                {
                    result[0] = WorkPrec_(0);
                    WorkPrec_ temp = manning_n_squared * pow(h, WorkPrec_(-7/3)) * sqrt(q1 * q1 + q2 * q2) * (-1);

                    result[1] = (temp * q1 - h * slope_x) * WorkPrec_(9.81);
                    result[2] = (temp * q2 - h * slope_y) * WorkPrec_(9.81);

                    return result;
                }
                else
                {
                    result[0] = WorkPrec_(0);
                    WorkPrec_ temp = manning_n_squared * pow(std::numeric_limits<WorkPrec_>::epsilon(), WorkPrec_(-7/3)) * sqrt(q1 * q1 + q2 * q2) * (-1);

                    result[1] = (temp * q1 -  std::numeric_limits<WorkPrec_>::epsilon() * slope_x) * WorkPrec_(9.81);
                    result[2] = (temp * q2 -  std::numeric_limits<WorkPrec_>::epsilon() * slope_y) * WorkPrec_(9.81);

                    return result;
                }

            }

        public:
            /**
             * \brief Implementation of source processing.
             *
             * \param vector The target vector.
             * \param bottom_slopes_x The vector containing the slope values (x direction).
             * \param bottom_slopes_y The vector containing the slope values.(y direction)
             * \param manning_n_squared The manning constant squared.
             *
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> value(DenseVector<WorkPrec_> & vector, DenseVector<WorkPrec_> & bottom_slopes_x, DenseVector<WorkPrec_> & bottom_slopes_y, WorkPrec_ manning_n_squared)
            {
                CONTEXT("When processing RelaxSolver source, type: simple");

                typename DenseVector<WorkPrec_>::ElementIterator writeelementiterator(vector.begin_elements());
                typename DenseVector<WorkPrec_>::ConstElementIterator bottomslopesxiterator(bottom_slopes_x.begin_elements());
                typename DenseVector<WorkPrec_>::ConstElementIterator bottomslopesyiterator(bottom_slopes_y.begin_elements());
                DenseVector<WorkPrec_> temp((unsigned long)(3), WorkPrec_(0));
                WorkPrec_ height, velocity1;
                for (typename DenseVector<WorkPrec_>::ConstElementIterator readelementiterator(vector.begin_elements()), vector_end(vector.end_elements()); readelementiterator != vector_end; ++readelementiterator, ++writeelementiterator, ++bottomslopesxiterator, ++bottomslopesyiterator)
                {
                    // Read height from given vector
                    height = *readelementiterator;
                    ++readelementiterator;

                    // Read velocity in x-direction
                    velocity1 = *readelementiterator;
                    ++readelementiterator;

                    // Y-velocity does not have to be fetched explicitly, it can be accessed through *readelementiterator

                    // Compute source-term for previously fetched values
                    temp = _source(height, velocity1, *readelementiterator, *bottomslopesxiterator, *bottomslopesyiterator, manning_n_squared);

                    // Write computed values back to the given vector
                    *writeelementiterator = temp[0]; // Write height

                    ++writeelementiterator;
                    *writeelementiterator = temp[1]; // Write x-velocity

                    ++writeelementiterator;
                    *writeelementiterator = temp[2]; // Write y-velocity
                }
                return vector;
            }
    };

    template <>
    struct SourceProcessing<source_types::SIMPLE, tags::CPU::SSE>
    {
        public:
            static DenseVector<float> value(DenseVector<float> & vector, DenseVector<float> & bottom_slopes_x, DenseVector<float> & bottom_slopes_y, float manning_n_squared);
            static DenseVector<double> value(DenseVector<double> & vector, DenseVector<double> & bottom_slopes_x, DenseVector<double> & bottom_slopes_y, double manning_n_squared);

    };
}
#endif
