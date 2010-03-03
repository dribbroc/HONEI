/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
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

#include <honei/util/exception.hh>

namespace honei
{
    /**
     * MatrixError is the base class for all Matrix related exceptions.
     *
     * \ingroup grpmatrixexceptions
     */
    class MatrixError :
        public Exception
    {
        protected:
            MatrixError(const std::string & message) throw ();
    };

    /**
     * MatrixColumnsDoNotMatch is thrown when the number of columns of a Matrix argument
     * does not match the expected number of columns.
     *
     * \ingroup grpmatrixexceptions
     */
    class MatrixColumnsDoNotMatch :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param columns Number of columns of the Matrix that arose the problem.
             * \param expected_columns Number of Columns that was expected by the operation.
             */
            MatrixColumnsDoNotMatch(unsigned long columns, unsigned long expected_columns) throw ();
    };

    /**
     * MatrixRowsDoNotMatch is thrown when the number of rows of a Matrix argument
     * does not match the expected number of rows.
     *
     * \ingroup grpmatrixexceptions
     */
    class MatrixRowsDoNotMatch :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param rows Number of rows of the Matrix that arose the problem.
             * \param expected_rows Number of rows that was expected by the operation.
             */
            MatrixRowsDoNotMatch(unsigned long rows, unsigned long expected_rows) throw ();
    };

    /**
     * MatrixSizeDoesNotMatch is thrown when the size of a square Matrix argument
     * does not match the expected size.
     *
     * \ingroup grpmatrixexceptions
     */
    class MatrixSizeDoesNotMatch :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param size Size of the square Matrix that arose the problem.
             * \param expected_size Size that was expected by the operation.
             */
            MatrixSizeDoesNotMatch(unsigned long size, unsigned long expected_size) throw ();
    };

    /**
     * MatrixIsNotSquare is thrown when a Matrix argument is not a square matrix, i.e. its
     * number of rows does not match its number of columns.
     *
     * \ingroup grpmatrixexceptions
     */
    class MatrixIsNotSquare :
        public MatrixError
    {
        public:
            /**
             * Constructor.
             *
             * \param rows Number of rows of the Matrix that arose the problem.
             * \param expected_columns Number fo columsn that was expected by the operation.
             */
            MatrixIsNotSquare(unsigned long rows, unsigned long expected_columns) throw ();
    };
}

#endif
