/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Sven Mallach <sven.mallach@honei.org>
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

#include <honei/cell/cell.hh>
#include <honei/libla/element_product.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float> &
    ElementProduct<tags::Cell>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When multiplying DenseMatrix<float> to DenseMatrix<float> (Cell):");

        if (a.rows() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        if (a.columns() != b.columns())
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());

        SPEFrameworkInstruction<2, float, rtm_dma> instruction(oc_element_product_dense_dense_float, a.elements(), b.elements(), a.rows() * a.columns());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        Matrix<float>::ConstElementIterator k(b.element_at(instruction.transfer_end()));
        for (MutableMatrix<float>::ElementIterator i(a.element_at(instruction.transfer_end())),
            i_end(a.end_elements()) ; i != i_end ; ++i, ++k)
        {
            *i *= *k;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }



    DenseVectorContinuousBase<float> &
    ElementProduct<tags::Cell>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
    {
        CONTEXT("When multiplying DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (Cell):");

        if (b.size() != a.size())
            throw VectorSizeDoesNotMatch(b.size(), a.size());

        SPEFrameworkInstruction<2, float, rtm_dma> instruction(oc_element_product_dense_dense_float, a.elements(), b.elements(), a.size());

        if (instruction.use_spe())
        {
            SPEManager::instance()->dispatch(instruction);
        }

        // Calculate the first 4 - a_offset % 4 elements on PPU.
        Vector<float>::ConstElementIterator j(b.begin_elements());
        for (Vector<float>::ElementIterator i(a.begin_elements()),
            i_end(a.element_at(instruction.transfer_begin())) ; i != i_end ; ++i, ++j)
        {
            *i *= *j;
        }

        Vector<float>::ConstElementIterator k(b.element_at(instruction.transfer_end()));
        for (Vector<float>::ElementIterator i(a.element_at(instruction.transfer_end())),
            i_end(a.end_elements()) ; i != i_end ; ++i, ++k)
        {
            *i *= *k;
        }

        if (instruction.use_spe())
            instruction.wait();

        return a;
    }

}
