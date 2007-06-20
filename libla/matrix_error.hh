/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_ERROR_HH
#define LIBLA_GUARD_MATRIX_ERROR_HH 1

#include <libutil/exception.hh>

/**
 * \file
 *
 * Declaration of exception classes which are related to Matrix.
 *
 * \ingroup grpmatrixexceptions
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A MatrixError is the base class for all Matrix related exceptions.
     *
     * \ingroup grpmatrixexceptions
     **/
    class MatrixError :
        public Exception
    {
        protected:
            MatrixError(const std::string & message) throw ();
    };

    /**
     * A MatrixColumnsDoNotMatch is thrown when the number of columns of a Matrix argument
     * does not match the expected number of columns.
     *
     * \ingroup grpmatrixexceptions
     **/
    class MatrixColumnsDoNotMatch :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param columns Amount of columns of the Matrix that arose the problem.
             * \param expected_columns Amount of Columns that was expected by the operation.
             **/
            MatrixColumnsDoNotMatch(unsigned long columns, unsigned long expected_columns) throw ();
    };

     /**
     * A MatrixRowsDoNotMatch is thrown when the number of rows of a Matrix argument
     * does not match the expected number of rows.
     *
     * \ingroup grpmatrixexceptions
     **/
    class MatrixRowsDoNotMatch :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param rows Amount of rows of the Matrix that arose the problem.
             * \param expected_rows Amount of rows that was expected by the operation.
             **/
            MatrixRowsDoNotMatch(unsigned long rows, unsigned long expected_rows) throw ();
    };

     /**
     * A MatrixMultiplicationError is thrown when the number of colums of a Matrix argument
     * does not match the number of rows of the Matrix argument to be multiplied with.
     *
     * \ingroup grpmatrixexceptions
     **/
    class MatrixMultiplicationError :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param columns Amount of columns of the first referenced matrix.
             * \param rows Amount of rows of the second referenced matrix.
             **/
            MatrixMultiplicationError(unsigned long columns, unsigned long rows) throw ();
    };

}

#endif
