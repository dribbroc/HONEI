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

#include <cell/cell.hh>
#include <libla/product.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace honei
{
    DenseVector<float>
    Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseVector<float> & x)
    {
        CONTEXT("When calculating DenseMatrix<float>-DenseVector<float> product (Cell):");

        if (x.size() != a.columns())
            throw VectorSizeDoesNotMatch(x.size(), a.columns());

        DenseVector<float> result(a.rows());

        Operand oa = { a.elements() };
        Operand ob = { x.elements() };
        Operand oc = { result.elements() };
        Operand od;
        od.u = a.rows();
        SPEInstruction instruction(oc_dense_dense_float_matrix_vector_product, x.size(), oa, ob, oc, od);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();

        return result;
    }

    DenseMatrix<float>
    Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When calculating DenseMatrix<float>-DenseMatrix<float> product (Cell):");

        if (a.columns() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.columns());

        DenseMatrix<float> result(a.rows(), b.columns(), float(0));

        Operand oa = { a.elements() };
        Operand ob = { b.elements() };
        Operand oc = { result.elements() };
        Operand od;
        Operand oe;
        od.u = a.rows();
        oe.u = b.columns();

        SPEInstruction instruction(oc_dense_dense_float_matrix_product, a.columns(), oa, ob, oc, od, oe);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();

        return result;
    }
}

