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

#include <libla/matrix_error.hh>
#include <libutil/stringify.hh>

#include <string>

using namespace pg512;

MatrixError::MatrixError(const std::string & message) throw () :
    Exception(message)
{
}

MatrixColumnsDoNotMatch::MatrixColumnsDoNotMatch(unsigned long columns, unsigned long expected_columns) throw () :
    MatrixError("Matrix column count '" + stringify(columns) + "' does not match expected column count '"
            + stringify(expected_columns) + "'")
{
}

MatrixRowsDoNotMatch::MatrixRowsDoNotMatch(unsigned long rows, unsigned long expected_rows) throw () :
    MatrixError("Matrix row count '" + stringify(rows) + "' does not match expected row count '"
            + stringify(expected_rows) + "'")
{
}

MatrixRowsDoNotMatch::MatrixMultiplicationError(unsigned long columns, unsigned long rows) throw () :
    MatrixError("First Matrix column count '" + stringify(columns) + "' does not match second matrix row count '"
            + stringify(rows) + "'")
{
}
