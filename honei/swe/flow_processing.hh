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

#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/swe/directions.hh>
/**
 * \file
 * Implementation of flow processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    template <typename Direction_, typename Tag_ = tags::CPU>
    struct FlowProcessing
    {
    };

    /**
     * \brief Implementation of flow processing in x direction.
     *
     * \ingroup grplibswe
     *
     **/
    template <typename Tag_>
    struct FlowProcessing<directions::X, Tag_>
    {
        private:
            /**
             * \brief Implementation of elementwise flow processing in x direction.
             * \param h The height at current position.
             * \param q1 The product of h and velocity in x direction.
             * \param q2 The product of h and velocity in y direction.
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
                static inline DenseVector<WorkPrec_> _flow_x(WorkPrec_ h, WorkPrec_ q1, WorkPrec_ q2)
                {
                    CONTEXT("When processing RelaxSolver flow, direction: X");

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
            /**
             * \brief Implementation of flow processing in x direction.
             * \param vector The vector for which the flow is to be computed.
             *
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> value(DenseVector<WorkPrec_> & vector)
            {
                CONTEXT("When processing RelaxSolver flow, direction: X");

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

    template <>
    struct FlowProcessing<directions::X, tags::Cell>
    {
        public:
            static DenseVector<float> & value(DenseVector<float> & vector);
            static DenseVector<double> & value(DenseVector<double> & vector)
            {
                FlowProcessing<directions::X, tags::CPU>::value(vector);
                return vector;
            }
    };

    /**
     * \brief Implementation of flow processing in y direction.
     *
     * \ingroup grplibswe
     *
     **/
    template <typename Tag_>
    struct FlowProcessing<directions::Y, Tag_>
    {
        private:
            /**
             * \brief Implementation of elementwise flow processing in x direction.
             * \param h The height at current position.
             * \param q1 The product of h and velocity in x direction.
             * \param q2 The product of h and velocity in y direction.
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> _flow_y(WorkPrec_ h, WorkPrec_ q1, WorkPrec_ q2)
            {
                CONTEXT("When processing RelaxSolver flow, direction: Y");

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
            /**
             * \brief Implementation of flow processing in x direction.
             * \param vector The vector for which the flow is to be computed.
             *
             * \ingroup grplibswe
             *
             **/
            template <typename WorkPrec_>
            static inline DenseVector<WorkPrec_> value(DenseVector<WorkPrec_> & vector)
            {
                CONTEXT("When processing RelaxSolver flow, direction: Y");

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

    template <>
    struct FlowProcessing<directions::Y, tags::Cell>
    {
        public:
            static DenseVector<float> & value(DenseVector<float> & vector);
            static DenseVector<double> & value(DenseVector<double> & vector)
            {
                FlowProcessing<directions::Y, tags::CPU>::value(vector);
                return vector;
            }
    };

///---------SSE-----------

    template <>
    struct FlowProcessing<directions::X, tags::CPU::SSE>
    {
        public:
            static DenseVector<float> value(DenseVector<float> & vector);
            static DenseVector<double> value(DenseVector<double> & vector);
    };

    template <>
    struct FlowProcessing<directions::Y, tags::CPU::SSE>
    {
        public:
            static DenseVector<float> value(DenseVector<float> & vector);
            static DenseVector<double> value(DenseVector<double> & vector);
    };

    template <>
    struct FlowProcessing<directions::X, tags::CPU::MultiCore::SSE>
    {
        public:
            static DenseVector<float> value(DenseVector<float> & vector)
            {
                return FlowProcessing<directions::X, tags::CPU::SSE>::value(vector);
            }
            static DenseVector<double> value(DenseVector<double> & vector)
            {
                return FlowProcessing<directions::X, tags::CPU::SSE>::value(vector);
            }
    };

    template <>
    struct FlowProcessing<directions::Y, tags::CPU::MultiCore::SSE>
    {
        public:
            static DenseVector<float> value(DenseVector<float> & vector)
            {
                return FlowProcessing<directions::Y, tags::CPU::SSE>::value(vector);
            }
            static DenseVector<double> value(DenseVector<double> & vector)
            {
                return FlowProcessing<directions::Y, tags::CPU::SSE>::value(vector);
            }
    };

}
#endif
