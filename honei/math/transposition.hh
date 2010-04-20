/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI MATH C++ library. Math is free software;
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

#ifndef MATH_GUARD_TRANSPOSITION_HH
#define MATH_GUARD_TRANSPOSITION_HH 1

#include <honei/la/sparse_matrix.hh>

namespace honei
{
    template <typename Tag_>
    struct Transposition
    {
        public:
        template <typename DT_>
        static inline void value(SparseMatrix<DT_> & source, SparseMatrix<DT_> & target)
        {
            CONTEXT("When transposing sparse matrix: ");
            if(source.rows() != target.columkns() || source.columns() != target.rows())
                throw InternalError("Inner matrix dimensions mismatch!");

            for(typename SparseMatrix<DT_>::NonZeroConstElementIterator i(source.begin_non_zero_elements()) ; i != source.end_non_zero_elements() ; ++i)
            {
                target(i->column(), i->row()) = *i;
            }
        }
    };
}

#endif
