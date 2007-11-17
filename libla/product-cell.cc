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

namespace bm_dv_product
{
    struct SPETask
    {
        honei::SPEInstruction spe_instruction;
        float * band;
        unsigned long start;
        unsigned long end;
    };
}

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

        unsigned int spe_count = 3;


        typedef std::vector<bm_dv_product::SPETask> TaskList;
        TaskList task_list;
        DenseVector<float> * results[4] = { new DenseVector<float>(b.size(), float(0)),
                                            new DenseVector<float>(b.size(), float(0)),
                                            new DenseVector<float>(b.size(), float(0)),
                                            new DenseVector<float>(b.size(), float(0)) };

        unsigned int counter = 0;

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

            Operand oa = { vi->elements() + quad_start };
            Operand ob = { b.elements() };
            Operand oc = { results[counter]->elements() + quad_start };
            Operand od, oe, of, og;
            od.u = quad_start;
            oe.u = quad_end;
            of.s = (signed)op_offset;
            og.u = op_offset % 4;
            bm_dv_product::SPETask task;
            task.band = vi->elements();
            task.start = start;
            task.end = end;
            task.spe_instruction= SPEInstruction(oc_banded_dense_float_matrix_vector_product, b.size(), oa, ob, oc, od, oe, of, og);
            task_list.push_back(task);

            /*for (unsigned long index = quad_end ; index < end ; index++) 
            {
                results[counter]->elements()[index] += vi->elements()[index] * b.elements()[index + op_offset];
            }*/
            counter = (counter + 1) % spe_count;
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
            Operand oa = { vi->elements() + quad_start};
            Operand ob = { b.elements() };
            Operand oc = { results[counter]->elements() + quad_start};
            Operand od, oe, of, og;
            od.u = quad_start;
            oe.u = quad_end;
            of.s = -1 * (signed)op_offset;
            og.u = (4 - (op_offset % 4)) % 4;
            bm_dv_product::SPETask task;
            task.band = vi->elements();
            task.start = start;
            task.end = end;
            task.spe_instruction= SPEInstruction(oc_banded_dense_float_matrix_vector_product, b.size(), oa, ob, oc, od, oe, of, og);
            task_list.push_back(task);

            for (unsigned long index = start ; index < quad_start ; index++)
            {
                results[counter]->elements()[index] += vi->elements()[index] * b.elements()[index - op_offset];
            }
            /*for (unsigned long index = quad_end ; index < end ; index++)
            {
                results[counter]->elements()[index] += vi->elements()[index] * b.elements()[index - op_offset];
            }*/
            counter = (counter + 1) % spe_count;
        }

        // dispatch instructions to spe's
        for (unsigned long i = 0 ; i < task_list.size() - (task_list.size() % spe_count) ; i = i + spe_count)
        {
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
                SPEManager::instance()->dispatch(task_list.at(i + j).spe_instruction);
            }
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
                for (unsigned long index = task_list.at(i + j).spe_instruction.instruction().e.u ; index < task_list.at(i).end ; index++) 
                {
                    results[j]->elements()[index] += 
                        task_list.at(i + j).band[index] *
                        b.elements()[index + task_list.at(i + j).spe_instruction.instruction().f.s];
                }
            }
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
            task_list.at(i + j).spe_instruction.wait();
            }
        }

        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            SPEManager::instance()->dispatch(task_list.at(i).spe_instruction);
        }
        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            for (unsigned long index = task_list.at(i).spe_instruction.instruction().e.u ; index < task_list.at(i).end ; index++)
            {
                ((float *)(task_list.at(i).spe_instruction.instruction().c.ea))[index] += 
                    task_list.at(i).band[index] * 
                    b.elements()[index + task_list.at(i).spe_instruction.instruction().f.s];
            }
        }
        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            task_list.at(i).spe_instruction.wait();
        }


        // 0 + 1
        // 2 + 3
        // 0 + 2
        Operand orc, ord;
        // hardcode transfer buffer size for now.
        orc.u = b.size() / (1024 * 4);
        ord.u = b.size() % (1024 * 4);
        ord.u &= ~0xF;

        unsigned result_rest_index(orc.u * 4096 + ord.u);

        ord.u *= 4;

        bool result_use_spe(true);

        if (0 == ord.u)
        {
            if (orc.u > 0)
            {
                ord.u = 16 * 1024;
            }
            else
            {
                result_use_spe = false;
            }
        }
        else
        {
            ++orc.u;
        }

        Operand or0 = { results[0]->elements() };
        Operand or1 = { results[1]->elements() };
        Operand or2 = { results[2]->elements() };
        Operand or3 = { results[3]->elements() };
        SPEInstruction result_x_instruction(oc_dense_dense_float_sum, 16 * 1024, or0, or1, orc, ord);
        SPEInstruction result_y_instruction(oc_dense_dense_float_sum, 16 * 1024, or2, or3, orc, ord);

        if (result_use_spe)
        {
            SPEManager::instance()->dispatch(result_x_instruction);
            SPEManager::instance()->dispatch(result_y_instruction);
        }

        Vector<float>::ConstElementIterator j1(results[1]->element_at(result_rest_index));
        Vector<float>::ConstElementIterator j3(results[3]->element_at(result_rest_index));
        Vector<float>::ElementIterator i0(results[0]->element_at(result_rest_index)), i0_end(results[0]->end_elements());
        Vector<float>::ElementIterator i2(results[2]->element_at(result_rest_index)), i2_end(results[2]->end_elements());
        for ( ; i0 != i0_end ; ++i0, ++j1)
        {
            *i0 += *j1;
        }
        for ( ; i2 != i2_end ; ++i2, ++j3)
        {
            *i2 += *j3;
        }

        if (result_use_spe)
        {
            result_x_instruction.wait();
            result_y_instruction.wait();
        }


        SPEInstruction result_z_instruction(oc_dense_dense_float_sum, 16 * 1024, or0, or2, orc, ord);

        if (result_use_spe)
        {
            SPEManager::instance()->dispatch(result_z_instruction);
        }

        Vector<float>::ConstElementIterator j2(results[2]->element_at(result_rest_index));
        i0 = results[0]->element_at(result_rest_index); i0_end = results[0]->end_elements();
        for ( ; i0 != i0_end ; ++i0, ++j2)
        {
            *i0 += *j2;
        }

        if (result_use_spe)
        {
            result_z_instruction.wait();
        }

        return *results[0];
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

