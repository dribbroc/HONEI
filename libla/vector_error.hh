/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_VECTOR_ERROR_HH
#define LIBLA_GUARD_VECTOR_ERROR_HH 1

#include <libutil/exception.hh>

#include <string>

/**
 * \file
 *
 * Declaration of exception classes which are related to Vector.
 *
 * \ingroup grpvectorexceptions
 **/
namespace honei
{
    /**
     * A VectorError is the base class for all Vector related exceptions.
     *
     * \ingroup grpvectorexceptions
     **/
    class VectorError :
        public Exception
    {
        protected:
            VectorError(const std::string & message) throw ();
    };

    /**
     * A VectorSizeDoesNotMatch is thrown when the size of a Vector argument
     * does not match the expected size.
     *
     * \ingroup grpvectorexceptions
     **/
    class VectorSizeDoesNotMatch :
        public VectorError
    {
        public:
            /**
             * Constructor.
             *
             * \param size Size of the vector that arose the problem.
             * \param expected_size Size that was expected by the operation.
             **/
            VectorSizeDoesNotMatch(unsigned long size, unsigned long expected_size) throw ();
    };
}

#endif
