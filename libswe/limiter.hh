/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Joachim Messer <joachim.messer@t-online.de>
 * Copyright (c) 2007 Volker Jung <volker.m.jung@t-online.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
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


#ifndef LIBSWE_GUARD_LIMITER_HH
#define LIBSWE_GUARD_LIMITER_HH 1

#include <libutil/tags.hh>

/**
 * \file
 *
 * Implementations for Limiter-based classes.
 *
 * \ingroup grplibswe
 **/

namespace pg512
{

    /**
     * \brief Limiter is the abstract base class for all Limiter-like classes used.
     *
     * \ingroup grplibswe
     **/

    template<typename Tag_ = tags::CPU> class Limiter
        {
            protected:

            /// Returns the maximum of two values.
            template<typename DataType_> static const DataType_ max(const DataType_ & left, const DataType_ & right)
            {
                return ( (left < right) ? right : left);
            }

            /// Returns the minimum of two values.            
            template<typename DataType_> static const DataType_ min(const DataType_ & left, const DataType_ & right)
            {
                return ( (left < right) ? left : right);
            }

            public:

            /// Abstract function for evaluating the limiterfunction.
            template<typename DataType_> static const DataType_ value(const DataType_ & scal);

        };


    /**
     * \brief MM_Limiter provides the MinMod-Limiter
     *
     * \ingroup grplibswe
     **/

    template<typename Tag_ = tags::CPU> class MM_Limiter : public Limiter<>
    {
        public:

        /// Evaluates due to the MinMod-Limiter specification
        template<typename DataType_> static const DataType_ value(const DataType_ & scal)
        {
            return( max<DataType_>( (DataType_) 0, min<DataType_>( (DataType_) 1, scal) ) );
        }
    };

    /**
     * \brief SB_Limiter provides the SuperBee-Limiter
     *
     * \ingroup grplibswe
     **/

    template<typename Tag_ = tags::CPU> class SB_Limiter : public Limiter<>
    {
        public:

        
        template<typename DataType_> static const DataType_ value(const DataType_ & scal)
        {
            return( max<DataType_>( (DataType_) 0, max<DataType_>( min<DataType_>(2*scal,(DataType_) 1), min<DataType_>(scal,(DataType_) 2) ) ) );
        }
    };

    /**
     * \brief MC_Limiter provides the Monotonized-Central-Limiter
     *
     * \ingroup grplibswe
     **/

    template<typename Tag_ = tags::CPU> class MC_Limiter : public Limiter<>
    {
        public:

        /// Evaluates due to the Monotonized-Central-Limiter specification
        template<typename DataType_> static const DataType_ value(const DataType_ & scal)
        {
            return( max<DataType_>( (DataType_) 0, min<DataType_>( min<DataType_>( 2*scal, (DataType_) 2), (1+scal)/2) ) );
        }
    };

    /**
     * \brief VL_Limiter provides the VanLeer-Limiter
     *
     * \ingroup grplibswe
     **/

    template<typename Tag_ = tags::CPU> class VL_Limiter : public Limiter<>
    {
        public:

        /// Evaluates due to the VanLeer-Limiter specification
        template<typename DataType_> static const DataType_ value(const DataType_ & scal)
        {
            return( (scal <= (DataType_) 0 ) ? (DataType_) 0 : (2*scal/(1+scal)) );
        }
    };

}

#endif
