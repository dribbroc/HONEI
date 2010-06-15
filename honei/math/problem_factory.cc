/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.tu-dortmund.de>
 *
 * This file is part of the HONEI math C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/math/problem_factory.hh>

namespace honei
{
    template class FileFactory<float, SparseMatrixELL<float>, SparseMatrixELL<float> >;
    template<> std::string FileFactory<float, SparseMatrixELL<float>, SparseMatrixELL<float> >::_filebase("4711");

    template class FileFactory<double, SparseMatrixELL<double>, SparseMatrixELL<double> >;
    template<> std::string FileFactory<double, SparseMatrixELL<double>, SparseMatrixELL<double> >::_filebase("4711");
}

