/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <iostream>

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

    DenseVector<float>
    Product<tags::Cell>::value(const BandedMatrix<float> & a, const DenseVector<float> & b)
    {
        CONTEXT("When calculating BandedMatrix<float>-DenseVector<float> product (Cell):");

        if (b.size() != a.columns())
            throw VectorSizeDoesNotMatch(b.size(), a.columns());

        DenseVector<float> result(a.rows(), float(0));

        unsigned long middle_index(a.rows() - 1);
        unsigned long quad_end, end, quad_start, start, x_offset, op_offset;


        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        for (BandedMatrix<float>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
                vi != vi_end ; ++vi)
        {
            if (! vi.exists())
                continue;

            op_offset = vi.index() - middle_index;
            start = 0;
            quad_start = 0;
            end = vi->size() - op_offset; //Calculation of the element-index to stop in iteration!
            quad_end = end - (end % 4);

            Operand oa = { vi->elements() };
            Operand ob = { b.elements() };
            Operand oc = { result.elements() };
            Operand od, oe, of, og;
            od.u = quad_start;
            oe.u = quad_end;
            of.s = (signed)op_offset;
            //og.u = (4- (op_offset % 4)) % 4;
            og.u = op_offset % 4;
            SPEInstruction instruction(oc_banded_dense_float_matrix_vector_product, b.size(), oa, ob, oc, od, oe, of, og);

            SPEManager::instance()->dispatch(instruction);
            instruction.wait();

            for (unsigned long index = quad_end ; index < end ; index++) 
            {
                result.elements()[index] += vi->elements()[index] * b.elements()[index + op_offset];
            }
        }

        // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        for (BandedMatrix<float>::ConstVectorIterator vi(a.begin_bands()), vi_end(a.band_at(middle_index)) ;
                vi != vi_end ; ++vi)
        {
            if (! vi.exists())
                continue;

            op_offset = middle_index - vi.index();
            start = op_offset; //Calculation of the element-index to start in iteration!
            quad_start = start + ((4 - (start % 4)) % 4);
            end = vi->size();
            quad_end = end - (end % 4);
            Operand oa = { vi->elements() };
            Operand ob = { b.elements() };
            Operand oc = { result.elements() };
            Operand od, oe, of, og;
            od.u = quad_start;
            oe.u = quad_end;
            of.s = -1 * (signed)op_offset;
            og.u = (4 - (op_offset % 4)) % 4;
            SPEInstruction instruction(oc_banded_dense_float_matrix_vector_product, b.size(), oa, ob, oc, od, oe, of, og);

            SPEManager::instance()->dispatch(instruction);
            instruction.wait();
            std::cout << "start "<<start<<" quad_start "<<quad_start<<std::endl;
            for (unsigned long index = start ; index < quad_start ; index++)
            {
                result.elements()[index] += vi->elements()[index] * b.elements()[index - op_offset];
            }
            for (unsigned long index = quad_end ; index < end ; index++)
            {
                result.elements()[index] += vi->elements()[index] * b.elements()[index - op_offset];
            }
        }

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

