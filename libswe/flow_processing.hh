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


#ifndef LIBSWE_GUARD_FLOW_PROCESSING_HH
#define LIBSWE_GUARD_FLOW_PROCESSING_HH 1

#include <libla/dense_vector.hh>
#include <libla/product.hh>
#include <libswe/directions.hh>
/**
 * \file
 * Implementation of flow processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    template <typename Tag_, typename Direction_>
    struct FlowProcessing
    {
    };

    /**
     * \brief Implementation of flow processing in x direction.
     *
     * \ingroup grplibswe
     *
     **/
    template <>
    struct FlowProcessing<tags::CPU, directions::X>
    {
        private:
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> _flow_x(WorkPrec_ h, WorkPrec_ q1, WorkPrec_ q2)
            {
                DenseVector<WorkPrec_> result((unsigned long)(3), WorkPrec_(0));

                if (fabs(h) >= std::numeric_limits<WorkPrec_>::epsilon())
                {
                    result[0] = q1;
                    result[1] = (q1 * q1 / h) + (WorkPrec_(9.81) * h * h / WorkPrec_(2));
                    result[2] = q1 * q2 / h;
                }
                else
                {
                    result[0] = q1;
                    result[1] = (q1 * q1 / std::numeric_limits<WorkPrec_>::epsilon());
                    result[2] = q1 * q2 / std::numeric_limits<WorkPrec_>::epsilon();
                }
                return result;

            }

        public:
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> value(DenseVector<WorkPrec_> & vector)
            {
                typename DenseVector<WorkPrec_>::ElementIterator writeelementiterator(vector.begin_elements());
                WorkPrec_ height, velocity1;
                DenseVector<WorkPrec_> temp((unsigned long)(3), WorkPrec_(0));
                for (typename DenseVector<WorkPrec_>::ConstElementIterator readelementiterator(vector.begin_elements()), vector_end(vector.end_elements()); readelementiterator != vector_end; ++readelementiterator, ++writeelementiterator)
                {
                    // Read height from given vector
                    height = *readelementiterator;
                    ++readelementiterator;

                    // Read x-velocity from given vector
                    velocity1 = *readelementiterator;
                    ++readelementiterator;

                    // Y-velocity does not have to be fetched explicitly, it can be accessed through *readelementiterator

                    // Compute flow on previously read values
                    temp = _flow_x<WorkPrec_>(height, velocity1, *readelementiterator);

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

    /**
     * \brief Implementation of flow processing in y direction.
     *
     * \ingroup grplibswe
     *
     **/
    template <>
    struct FlowProcessing<tags::CPU, directions::Y>
    {
        private:
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> _flow_y(WorkPrec_ h, WorkPrec_ q1, WorkPrec_ q2)
            {
                DenseVector<WorkPrec_> result((unsigned long)(3), WorkPrec_(0));

                if (fabs(h) >= std::numeric_limits<WorkPrec_>::epsilon())
                {
                    result[0] = q2;
                    result[1] = q1 * q2 / h;
                    result[2] = (q2 * q2 / h) + (WorkPrec_(9.81) * h * h / WorkPrec_(2));
                }
                else
                {
                    result[0] = q2;
                    result[1] = q1 * q2 / std::numeric_limits<WorkPrec_>::epsilon();
                    result[2] = q2 * q2 / std::numeric_limits<WorkPrec_>::epsilon();
                }
                return result;

            }

        public:
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> value(DenseVector<WorkPrec_> & vector)
            {
                typename DenseVector<WorkPrec_>::ElementIterator writeelementiterator(vector.begin_elements());
                WorkPrec_ height, velocity1;
                DenseVector<WorkPrec_> temp((unsigned long)(3), WorkPrec_(0));
                for (typename DenseVector<WorkPrec_>::ConstElementIterator readelementiterator(vector.begin_elements()), vector_end(vector.end_elements()); readelementiterator != vector_end; ++readelementiterator, ++writeelementiterator)
                {
                    // Read height from given vector
                    height = *readelementiterator;
                    ++readelementiterator;
                    // Read x-velocity from given vector
                    velocity1 = *readelementiterator;
                    ++readelementiterator;

                    // Y-velocity does not have to be fetched explicitly, it can be accessed through *readelementiterator

                    // Compute flow on previously read values
                    temp = _flow_y<WorkPrec_>(height, velocity1, *readelementiterator);

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

    }
#endif