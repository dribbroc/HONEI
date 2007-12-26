/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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
#include <libla/sum.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace bm_dv_product
{
    struct SPETask
    {
        honei::SPEInstruction spe_instruction;
        float * band;
        unsigned int result_counter;
        unsigned long start;
        unsigned long end;
        unsigned long quad_start;
        unsigned long quad_end;
        bool use_spe;
    };
}

namespace honei
{
    using namespace cell;

    DenseVector<float> Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseVector<float> & x)
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

    DenseVector<float> Product<tags::Cell>::value(const BandedMatrix<float> & a, const DenseVector<float> & b)
    {
        CONTEXT("When calculating BandedMatrix<float>-DenseVector<float> product (Cell):");

        if (b.size() != a.columns())
            throw VectorSizeDoesNotMatch(b.size(), a.columns());

        unsigned int spe_count (std::min((unsigned)4 , SPEManager::instance()->spe_count()));
        //unsigned int spe_count(1);

        typedef std::vector<bm_dv_product::SPETask> TaskList;
        TaskList task_list;
        DenseVector<float> * results[4] = {
            new DenseVector<float>(b.size(), float(0)),
            new DenseVector<float>(b.size(), float(0)),
            new DenseVector<float>(b.size(), float(0)),
            new DenseVector<float>(b.size(), float(0)) };

        unsigned int counter = 0;

        unsigned long middle_index(a.rows() - 1);
        unsigned long quad_end, end, quad_start, start, x_offset, op_offset;


        for (BandedMatrix<float>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
                band != band_end ; ++band)
        {
            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            if (band.index() >= middle_index)
            {

                op_offset = band.index() - middle_index;
                start = 0;
                quad_start = 0;
                end = band->size() - op_offset; //Calculation of the element-index to stop in iteration!
                quad_end = end - (end % 4);

                Operand oa = { band->elements() + quad_start };
                Operand ob = { b.elements() + quad_start + op_offset - (op_offset % 4) };
                Operand oc = { results[counter]->elements() + quad_start };
                Operand od, oe, of, og;

                /// \todo use such a transfer size, that od.u is at least 2
                od.u = (quad_end - quad_start) / (1000 * 4);
                oe.u = (quad_end - quad_start) % (1000 * 4);
                if (0 == oe.u)
                {
                    if (od.u > 0)
                    {
                        oe.u = 1000 * 4;
                    }
                }
                else
                {
                    ++od.u;
                }

                /// \todo Find a way to use vector size above signed long max size
                of.s = (signed)op_offset;
                og.u = op_offset % 4;
                bm_dv_product::SPETask task;
                task.result_counter = counter;
                task.band = band->elements();
                task.start = start;
                task.end = end;
                if(quad_end > quad_start)
                {
                    task.use_spe = true;
                }
                else
                {
                    task.use_spe = false;
                    quad_start = 0;
                    quad_end = 0;
                }
                task.quad_start = quad_start;
                task.quad_end = quad_end;
                //std::cout<< "above D: "<<od.u<<" E: "<<oe.u <<" spe: "<<task.use_spe<<" q_e: "<<quad_end<<std::endl;
                task.spe_instruction = SPEInstruction(oc_banded_dense_float_matrix_vector_product, 1000 * 4, oa, ob, oc, od, oe, of, og);
                task_list.push_back(task);

                /*for (unsigned long index = quad_end ; index < end ; index++) 
                  {
                  results[counter]->elements()[index] += vi->elements()[index] * b.elements()[index + op_offset];
                  }*/
                counter = (counter + 1) % spe_count;
            }

            // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
            else
            {

                op_offset = middle_index - band.index();
                start = op_offset; //Calculation of the element-index to start in iteration!
                quad_start = start + ((4 - (start % 4)) % 4);
                end = band->size();
                quad_end = end - (end % 4);
                Operand oa = { band->elements() + quad_start};
                Operand ob = { b.elements() + quad_start - op_offset - ((4 - (op_offset % 4)) % 4)};
                Operand oc = { results[counter]->elements() + quad_start};
                Operand od, oe, of, og;

                od.u = (quad_end - quad_start) / (1000 * 4);
                oe.u = (quad_end - quad_start) % (1000 * 4);
                if (0 == oe.u)
                {
                    if (od.u > 0)
                    {
                        oe.u = 1000 * 4;
                    }
                }
                else
                {
                    ++od.u;
                }

                of.s = -1 * (signed)op_offset;
                og.u = (4 - (op_offset % 4)) % 4;
                bm_dv_product::SPETask task;
                task.result_counter = counter;
                task.band = band->elements();
                task.start = start;
                task.end = end;
                if(quad_end > quad_start)
                {
                    task.use_spe = true;
                }
                else
                {
                    task.use_spe = false;
                    quad_start = 0;
                    quad_end = start;
                }
                task.quad_start = quad_start;
                task.quad_end = quad_end;
                //std::cout<< "below D: "<<od.u<<" E: "<<oe.u <<" spe: "<<task.use_spe<<" q_e: "<<quad_end<<std::endl;
                task.spe_instruction = SPEInstruction(oc_banded_dense_float_matrix_vector_product, 1000 * 4, oa, ob, oc, od, oe, of, og);
                task_list.push_back(task);

                for (unsigned long index = start ; index < quad_start ; index++)
                {
                    results[counter]->elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                }
                /*for (unsigned long index = quad_end ; index < end ; index++)
                  {
                  results[counter]->elements()[index] += vi->elements()[index] * b.elements()[index - op_offset];
                  }*/
                counter = (counter + 1) % spe_count;
            }
        }

        // dispatch instructions to spe's
        for (unsigned long i = 0 ; i < task_list.size() - (task_list.size() % spe_count) ; i = i + spe_count)
        {
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
                if(task_list.at(i + j).use_spe) SPEManager::instance()->dispatch(task_list.at(i + j).spe_instruction);
            }
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
                for (unsigned long index = task_list.at(i + j).quad_end ; index < task_list.at(i + j).end ; index++) 
                {
                    results[task_list.at(i + j).result_counter]->elements()[index] +=
                        task_list.at(i + j).band[index] *
                        b.elements()[index + task_list.at(i + j).spe_instruction.instruction().f.s];
                }
            }
            for (unsigned long j = 0 ; j < spe_count ; j++)
            {
                if(task_list.at(i + j).use_spe) task_list.at(i + j).spe_instruction.wait();
            }
        }

        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            if(task_list.at(i).use_spe) SPEManager::instance()->dispatch(task_list.at(i).spe_instruction);
        }
        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            for (unsigned long index = task_list.at(i).quad_end ; index < task_list.at(i).end ; index++)
            {
                results[task_list.at(i).result_counter]->elements()[index] +=
                    task_list.at(i).band[index] *
                    b.elements()[index + task_list.at(i).spe_instruction.instruction().f.s];
            }
        }
        for (unsigned long i = task_list.size() - (task_list.size() % spe_count) ; i < task_list.size() ; ++i)
        {
            if(task_list.at(i).use_spe) task_list.at(i).spe_instruction.wait();
        }


        // 0 + 1
        // 2 + 3
        // 0 + 2
        if (spe_count > 1)
        {
            Sum<tags::Cell>::value(*results[0], *results[1]);
            Sum<tags::Cell>::value(*results[2], *results[3]);
            Sum<tags::Cell>::value(*results[0], *results[2]);
        }
        DenseVector<float> result (*results[0]);
        for (unsigned long i(0) ; i < 4 ; ++i)
        {
            delete results[i];
        }
        return result;
    }

    DenseMatrix<float> Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When calculating DenseMatrix<float>-DenseMatrix<float> product (Cell):");

        if (a.columns() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.columns());

        /* \\\todo: Tiling. */

        /* We assume to have complete matrices or complete matrix tiles
         * at this point. So later here will be another for-loop that calls the
         * spe program for every pair (triple with result) of tiles.
         * This for loop will be responsible for transfering the logically
         * correct pairs of tiles to the spe.
         * The spe program itself can then compute the same way, as if it would
         * do, if it got complete matrices as operands.
         */

        unsigned a_tile_rows;
        unsigned a_tile_columns;
        unsigned b_tile_rows;
        unsigned b_tile_columns;
        unsigned r_tile_rows;
        unsigned r_tile_columns;

        // Beginning of double buffering.

        /* For the first source Matrix (a) and the result Matrix (r) we apply double buffering,
         * while we always transfer the whole second source Matrix (b).
         * For a and r we always pull as many complete rows that fit into 16 KB.
         */

        a_tile_rows = a.rows(); // Later use tiles' rows here
        a_tile_columns = a.columns(); // Later use tiles' columns here
        b_tile_rows = b.rows(); // Later user tiles' rows here
        b_tile_columns = b.columns(); // Later use tiles' columns here

        DenseMatrix<float> result(a_tile_rows, b_tile_columns, 0.0f);

        r_tile_rows = result.rows(); // Later use tiles' rows here
        r_tile_columns = result.columns(); // Later use tiles' columns here.

        Operand oa = { a.elements() }; // Later use tiles' elements here
        Operand ob = { b.elements() }; // Later use tiles' elements here
        Operand oc = { result.elements() }; // Later use tiles' elements here
        Operand od, oe;
        od.u = a.columns(); // Later use tiles' columns here
        oe.u = b.columns(); // Later use tiles' columns here

        // We now need to know how much of the tiles we can transfer at most
        // using one dma transfer.
        // For the tile of a the number of columns decides, because we can then see
        // how many rows we can transfer within 16 KB. We do the same for the tile of r.
        // There is only one problem. Because r can have less or more columns as a, we
        // can have different default transfer sizes for them. But we need exactly the
        // same number of rows of a and r to let the spe-program work correctly.
        // So we need a # of rows, for that both transfersizes are % 16 == 0.

        Operand of; // Will transfer the number of full transfers for a.
        Operand og; // Will transfer the number of blocks to reserve for b.
        Operand oh; // Will transfer the size of a last transfer for a,
        Operand oi; // Will transfer the size of a last transfer for b.
        Operand oj; // Will transfer the size of a last transfer for r.
        Operand ok; // Will transfer the default transfer size for r.

        unsigned a_row_bytes(a_tile_columns * 4); // Bytes per row
        unsigned a_rows_per_transfer(16384 / a_row_bytes);

        unsigned r_row_bytes(r_tile_columns * 4);
        unsigned r_rows_per_transfer(16384 / r_row_bytes);

        // r can have more or less columns than a.
        // We must start with the bigger one to assure that we don't
        // exceed 16 KB.
        unsigned rows_per_transfer;
        a_tile_columns > r_tile_columns ? rows_per_transfer = a_rows_per_transfer : rows_per_transfer = r_rows_per_transfer;
        unsigned r_def_t_size(rows_per_transfer * r_row_bytes);
        unsigned a_def_t_size(rows_per_transfer * a_row_bytes);

        for ( ; ; ) // Need two default aligned transfer size as closest to 16384 as possible.
        {
            if (a_def_t_size % 16 == 0 && r_def_t_size % 16 == 0)
            {
                break;
            }
            a_def_t_size -= a_row_bytes;
            r_def_t_size -= r_row_bytes;
            rows_per_transfer--;
            if (a_def_t_size <= 0 || r_def_t_size <= 0)
            {
                std::cout << "Could not find an alignment" << std::endl;
            }
        }

        unsigned b_bytes(b_tile_rows * b_tile_columns * 4);
        of.u = a_tile_rows / rows_per_transfer;
        og.u = b_bytes / 16384;
        oh.u = (a_tile_rows % rows_per_transfer) * a_row_bytes;
        oi.u = b_bytes % 16384;
        oj.u = (r_tile_rows % rows_per_transfer) * r_row_bytes;
        ok.u = r_def_t_size;

        // If there is a rest for a and r to transfer, we need one more transfer...
        // Else if there is no rest, the last transfer size is equal to the normal one.
        if (0 != oh.u || 0 != oj.u)
        {
            of.u++;
        }
        else
        {
            oh.u = a_def_t_size;
            oj.u = r_def_t_size;
        }

        // If a fits into one transfer, we should say that it is one transfer needed and set last size.
        if (0 == of.u)
        {
            of.u++;
            oh.u = a_row_bytes * a_tile_columns;
        }

        if (oi.u != 0) // if we need one more block for rest rows of b.
        {
            og.u++;
        }

        SPEInstruction instruction(oc_dense_dense_float_matrix_product, a_def_t_size, oa, ob, oc, od, oe, of, og, oh, oi, oj, ok);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();

        return result;
    }

    DenseMatrix<float>
        Product<tags::Cell>::value(const SparseMatrix<float> & a, const DenseMatrix<float> & b)
        {
            CONTEXT("When calculating SparseMatrix<float>-DenseMatrix<float> product (Cell):");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<float> result(a.rows(), b.columns(), 0.0f);

            for(SparseMatrix<float>::ConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                ScaledSum<tags::Cell>::value(result[i.row()], b[i.column()], *i);
            }

            return result;
        }
}
