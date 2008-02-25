/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/libla/product.hh>
#include <honei/libla/sum.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/profiler.hh>
#include <honei/libutil/spe_transfer_list.hh>
#include <list>

//#include <honei/libutil/time_stamp.hh>

namespace honei
{
    using namespace cell;

    DenseVector<float> Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseVector<float> & x)
    {
        CONTEXT("When calculating DenseMatrix<float>-DenseVector<float> product (Cell):");

        if (x.size() != a.columns())
            throw VectorSizeDoesNotMatch(x.size(), a.columns());

        DenseVector<float> result(a.rows(), 0.0f);
        std::list<SPEInstruction *> instructions;

        unsigned long spe_count(Configuration::instance()->get_value("cell::product_dense_matrix_dense_vector_float", 4ul));
        spe_count = std::min(spe_count, SPEManager::instance()->spe_count());
        unsigned rows_per_spe(a.rows() / spe_count);
        unsigned ppu_rows(0);
        bool use_spe(true);

        if (a.rows() < (4 * spe_count))
        {
            ppu_rows = a.rows();
            spe_count = 0;
            rows_per_spe = 0;
            use_spe = false;
        }

        if (use_spe)
        {
            unsigned parts(spe_count);
            while (rows_per_spe > 4096)
            {
                rows_per_spe /= 2;
                parts *= 2;
            }
            rows_per_spe & ~0x3;
            ppu_rows = a.rows() - (spe_count * rows_per_spe);

            union addr
            {
                float * ptr;
                unsigned long long value;
            };

            unsigned a_t_size(a.columns() * (4096 / (a.columns())));
            while (a_t_size % 4 != 0)
            {
                a_t_size -= a.columns();
            }
            a_t_size *= 4;

            for (unsigned i(0) ; i < parts ; i++)
            {
                addr one = { (a.elements() +  a.columns()) };
                unsigned a_offset(one.value & 0xf);

                Operand oa = { a.elements() + (rows_per_spe * i * a.columns()) };
                Operand ob = { x.elements() };
                Operand oc = { result.elements() + (rows_per_spe * i) };

                Operand od, oe, of, og;
                // Transfer only full rows.
                od.u = (4 * rows_per_spe * a.columns()) / a_t_size;
                oe.u = (4 * rows_per_spe * a.columns()) % a_t_size;
                of.u = (4 * x.size()) / 16384;
                og.u = (4 * x.size()) % 16384;

                if (0 == oe.u)
                {
                    if (od.u > 0)
                    {
                        oe.u = a_t_size;
                    }
                }
                else
                {
                    ++od.u;
                }

                if (0 == og.u)
                {
                    if (of.u > 0)
                    {
                        og.u = 16384;
                    }
                }
                else
                {
                    ++of.u;
                }

                Operand oh;
                oh.u = a.columns();
                Operand oi;
                oi.u = a_offset / sizeof(float);

                if (i % spe_count == 0) // wait for the last parts to have finished.
                {
                    for(std::list<SPEInstruction *>::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; i++)
                    {
                        (*i)->wait();
                        delete *i;
                        instructions.erase(i);
                    }
                }

                SPEInstruction * instruction = new SPEInstruction(oc_product_dense_matrix_dense_vector_float, a_t_size, oa, ob, oc, od, oe, of, og, oh, oi);
                instructions.push_back(instruction);
                SPEManager::instance()->dispatch(*instruction);
            }
        }

        for (unsigned i(0) ; i < ppu_rows ; i++)
        {
            DenseVectorRange<float> row = a[(spe_count * rows_per_spe) + i];
            for (Vector<float>::ConstElementIterator c(row.begin_elements()), c_end(row.end_elements()), d(x.begin_elements()) ; c != c_end ; ++c, ++d)
            {
                result[(spe_count * rows_per_spe) + i] += *c * *d;
            }
        }

        if (use_spe)
        {
            for(std::list<SPEInstruction *>::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; i++)
            {
                (*i)->wait();
                delete *i;
            }
        }

        return result;
    }

    DenseVector<float> Product<tags::Cell>::value(const BandedMatrix<float> & a, const DenseVector<float> & b)
    {
        CONTEXT("When calculating BandedMatrix<float>-DenseVector<float> product (Cell):");

        PROFILER_START("Product<Cell>::value(bm, dv)");
        if (b.size() != a.columns())
            throw VectorSizeDoesNotMatch(b.size(), a.columns());

        /// \todo Remove SPEIQ list when one SPEInstructionQueue can handle more than 8 instructions.
        std::list<SPEInstructionQueue *> iq_upper, iq_lower;
        //TimeStamp dt,ct, as1, as2;
        //ct.take();
        PROFILER_START("Product<Cell>::value(bm, dv)->DV(size, 0)");
        DenseVector<float> result(b.size(), 0.0f);
        PROFILER_STOP("Product<Cell>::value(bm, dv)->DV(size, 0)");
        //dt.take();
        //std::cout<<"dv(0): "<<dt.sec() - ct.sec() << " "<<dt.usec() - ct.usec()<<std::endl;
        /// \todo Fill the result vector on the spu side.
        //DenseVector<float> result(b.size());
        //as1.take();

        unsigned long middle_index(a.rows() - 1);
        unsigned long quad_end, end, quad_start, start, x_offset, op_offset;
        for (BandedMatrix<float>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
                band != band_end ; ++band)
        {
            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            if (band.index() >= middle_index)
            {
                {
                    // Lower result part
                    op_offset = band.index() - middle_index;
                    start = 0;
                    quad_start = 0;
                    end = std::min(band->size() - op_offset, result.size() / 2); //Calculation of the element-index to stop in iteration!
                    quad_end = end - (end % 4);

                    Operand oa = { band->elements() + quad_start };
                    Operand ob = { b.elements() + quad_start + op_offset - (op_offset % 4) };
                    Operand oc = { result.elements() + quad_start };
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

                    og.u = op_offset % 4;
                    if(quad_end > quad_start)
                    {
                        if (iq_lower.empty() || iq_lower.back()->size() == 7)
                        {
                            iq_lower.push_back(new SPEInstructionQueue);
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = 0;
                    }

                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                }

                {
                    // Upper result part
                    op_offset = band.index() - middle_index;
                    start = result.size() / 2;
                    quad_start = start + ((4 - (start % 4)) % 4);
                    end = band->size() - op_offset; //Calculation of the element-index to stop in iteration!
                    quad_end = end - (end % 4);

                    Operand oa = { band->elements() + quad_start };
                    Operand ob = { b.elements() + quad_start + op_offset - (op_offset % 4) };
                    Operand oc = { result.elements() + quad_start };
                    Operand od, oe, of, og, oh;

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

                    og.u = op_offset % 4;
                    if(quad_end > quad_start)
                    {
                        if (iq_upper.empty() || iq_upper.back()->size() == 7)
                        {
                            iq_upper.push_back(new SPEInstructionQueue);
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = result.size() / 2;
                        quad_end = result.size() / 2;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                }
            }

            // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
            else
            {
                {
                    // Lower result part
                    op_offset = middle_index - band.index();
                    start = op_offset; //Calculation of the element-index to start in iteration!
                    quad_start = start + ((4 - (start % 4)) % 4);
                    end = result.size() / 2;
                    quad_end = end - (end % 4);
                    Operand oa = { band->elements() + quad_start};
                    Operand ob = { b.elements() + quad_start - op_offset - ((4 - (op_offset % 4)) % 4)};
                    Operand oc = { result.elements() + quad_start};
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

                    og.u = (4 - (op_offset % 4)) % 4;
                    if(quad_end > quad_start)
                    {
                        if (iq_lower.empty() || iq_lower.back()->size() == 7)
                        {
                            iq_lower.push_back(new SPEInstructionQueue);
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = start;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                }

                {
                    // Upper result part
                    op_offset = middle_index - band.index();
                    start = std::max(op_offset, result.size() / 2); //Calculation of the element-index to start in iteration!
                    quad_start = start + ((4 - (start % 4)) % 4);
                    end = band->size();
                    quad_end = end - (end % 4);
                    Operand oa = { band->elements() + quad_start};
                    Operand ob = { b.elements() + quad_start - op_offset - ((4 - (op_offset % 4)) % 4)};
                    Operand oc = { result.elements() + quad_start};
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

                    og.u = (4 - (op_offset % 4)) % 4;
                    if(quad_end > quad_start)
                    {
                        if (iq_upper.empty() || iq_upper.back()->size() == 7)
                        {
                            iq_upper.push_back(new SPEInstructionQueue);
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_float, 1000 * 4, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = start;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                }
            }
        }
        /*as2.take();
        std::cout<<"assembly: "<<as2.sec() - as1.sec() << " "<<as2.usec() - as1.usec()<<std::endl;
        TimeStamp at, bt;
        at.take();*/
        PROFILER_START("Product<Cell>::value(bm, dv)->dispatch");
        for (std::list<SPEInstructionQueue *>::iterator i(iq_lower.begin()), i_end(iq_lower.end()), j(iq_upper.begin()), j_end(iq_upper.end()) ;
                (i != i_end || j != j_end) ; )
        {
            if (i != i_end) SPEManager::instance()->dispatch(*(*i));
            if (j != j_end) SPEManager::instance()->dispatch(*(*j));
            /// \todo calc scalar parts here
            if (i != i_end) (*i)->wait();
            if (j != j_end) (*j)->wait();
            if (i != i_end) ++i;
            if (j != j_end) ++j;
        }
        for (std::list<SPEInstructionQueue *>::iterator i(iq_lower.begin()), i_end(iq_lower.end()), j(iq_upper.begin()), j_end(iq_upper.end()) ;
                (i != i_end || j != j_end) ; )
        {
            if (i != i_end) delete *i;
            if (j != j_end) delete *j;
            if (i != i_end) ++i;
            if (j != j_end) ++j;
        }
        //bt.take();
        //std::cout<<"wait: "<<bt.sec() - at.sec() << " "<<bt.usec() - at.usec()<<std::endl<<std::endl;
        PROFILER_STOP("Product<Cell>::value(bm, dv)->dispatch");
        PROFILER_STOP("Product<Cell>::value(bm, dv)");

        return result;
    }

    DenseVector<double> Product<tags::Cell>::value(const BandedMatrix<double> & a, const DenseVector<double> & b)
    {
        CONTEXT("When calculating BandedMatrix<double>-DenseVector<double> product (Cell):");

        if (b.size() != a.columns())
            throw VectorSizeDoesNotMatch(b.size(), a.columns());

        /// \todo Remove SPEIQ list when one SPEInstructionQueue can handle more than 8 instructions.
        std::list<SPEInstructionQueue *> iq_upper, iq_lower;
        //TimeStamp dt,ct, as1, as2;
        //ct.take();
        DenseVector<double> result(b.size(), 0.0f);
        //dt.take();
        //std::cout<<"dv(0): "<<dt.sec() - ct.sec() << " "<<dt.usec() - ct.usec()<<std::endl;
        /// \todo Fill the result vector on the spu side.
        //DenseVector<double> result(b.size());
        //as1.take();

        unsigned long middle_index(a.rows() - 1);
        unsigned long quad_end, end, quad_start, start, x_offset, op_offset;
        for (BandedMatrix<double>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
                band != band_end ; ++band)
        {
            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            if (band.index() >= middle_index)
            {
                {
                    // Lower result part
                    op_offset = band.index() - middle_index;
                    start = 0;
                    quad_start = 0;
                    end = std::min(band->size() - op_offset, result.size() / 2); //Calculation of the element-index to stop in iteration!
                    quad_end = end - (end % 2);

                    Operand oa = { band->elements() + quad_start };
                    Operand ob = { b.elements() + quad_start + op_offset - (op_offset % 2) };
                    Operand oc = { result.elements() + quad_start };
                    Operand od, oe, of, og;
                    /// \todo use such a transfer size, that od.u is at least 2
                    od.u = (quad_end - quad_start) / (1000 * 2);
                    oe.u = (quad_end - quad_start) % (1000 * 2);
                    if (0 == oe.u)
                    {
                        if (od.u > 0)
                        {
                            oe.u = 1000 * 2;
                        }
                    }
                    else
                    {
                        ++od.u;
                    }

                    og.u = op_offset % 2;
                    if(quad_end > quad_start)
                    {
                        if (iq_lower.empty() || iq_lower.back()->size() == 7)
                        {
                            iq_lower.push_back(new SPEInstructionQueue);
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = 0;
                    }

                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                }

                {
                    // Upper result part
                    op_offset = band.index() - middle_index;
                    start = result.size() / 2;
                    quad_start = start + ((2 - (start % 2)) % 2);
                    end = band->size() - op_offset; //Calculation of the element-index to stop in iteration!
                    quad_end = end - (end % 2);

                    Operand oa = { band->elements() + quad_start };
                    Operand ob = { b.elements() + quad_start + op_offset - (op_offset % 2) };
                    Operand oc = { result.elements() + quad_start };
                    Operand od, oe, of, og, oh;

                    /// \todo use such a transfer size, that od.u is at least 2
                    od.u = (quad_end - quad_start) / (1000 * 2);
                    oe.u = (quad_end - quad_start) % (1000 * 2);
                    if (0 == oe.u)
                    {
                        if (od.u > 0)
                        {
                            oe.u = 1000 * 2;
                        }
                    }
                    else
                    {
                        ++od.u;
                    }

                    og.u = op_offset % 2;
                    if(quad_end > quad_start)
                    {
                        if (iq_upper.empty() || iq_upper.back()->size() == 7)
                        {
                            iq_upper.push_back(new SPEInstructionQueue);
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = result.size() / 2;
                        quad_end = result.size() / 2;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index + op_offset];
                    }
                }
            }

            // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
            else
            {
                {
                    // Lower result part
                    op_offset = middle_index - band.index();
                    start = op_offset; //Calculation of the element-index to start in iteration!
                    quad_start = start + ((2 - (start % 2)) % 2);
                    end = result.size() / 2;
                    quad_end = end - (end % 2);
                    Operand oa = { band->elements() + quad_start};
                    Operand ob = { b.elements() + quad_start - op_offset - ((2 - (op_offset % 2)) % 2)};
                    Operand oc = { result.elements() + quad_start};
                    Operand od, oe, of, og;

                    od.u = (quad_end - quad_start) / (1000 * 2);
                    oe.u = (quad_end - quad_start) % (1000 * 2);
                    if (0 == oe.u)
                    {
                        if (od.u > 0)
                        {
                            oe.u = 1000 * 2;
                        }
                    }
                    else
                    {
                        ++od.u;
                    }

                    og.u = (2 - (op_offset % 2)) % 2;
                    if(quad_end > quad_start)
                    {
                        if (iq_lower.empty() || iq_lower.back()->size() == 7)
                        {
                            iq_lower.push_back(new SPEInstructionQueue);
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_lower.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = start;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                }

                {
                    // Upper result part
                    op_offset = middle_index - band.index();
                    start = std::max(op_offset, result.size() / 2); //Calculation of the element-index to start in iteration!
                    quad_start = start + ((2 - (start % 2)) % 2);
                    end = band->size();
                    quad_end = end - (end % 2);
                    Operand oa = { band->elements() + quad_start};
                    Operand ob = { b.elements() + quad_start - op_offset - ((2 - (op_offset % 2)) % 2)};
                    Operand oc = { result.elements() + quad_start};
                    Operand od, oe, of, og;

                    od.u = (quad_end - quad_start) / (1000 * 2);
                    oe.u = (quad_end - quad_start) % (1000 * 2);
                    if (0 == oe.u)
                    {
                        if (od.u > 0)
                        {
                            oe.u = 1000 * 2;
                        }
                    }
                    else
                    {
                        ++od.u;
                    }

                    og.u = (2 - (op_offset % 2)) % 2;
                    if(quad_end > quad_start)
                    {
                        if (iq_upper.empty() || iq_upper.back()->size() == 7)
                        {
                            iq_upper.push_back(new SPEInstructionQueue);
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                        else
                        {
                            iq_upper.back()->push_back(SPEInstruction(oc_product_banded_matrix_dense_vector_double, 1000 * 2, oa, ob, oc, od, oe, of, og));
                        }
                    }
                    else
                    {
                        quad_start = 0;
                        quad_end = start;
                    }

                    for (unsigned long index = start ; index < quad_start ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                    for (unsigned long index = quad_end ; index < end ; index++)
                    {
                        result.elements()[index] += band->elements()[index] * b.elements()[index - op_offset];
                    }
                }
            }
        }
        /*as2.take();
        std::cout<<"assembly: "<<as2.sec() - as1.sec() << " "<<as2.usec() - as1.usec()<<std::endl;
        TimeStamp at, bt;
        at.take();*/
        for (std::list<SPEInstructionQueue *>::iterator i(iq_lower.begin()), i_end(iq_lower.end()), j(iq_upper.begin()), j_end(iq_upper.end()) ;
                (i != i_end || j != j_end) ; )
        {
            if (i != i_end) SPEManager::instance()->dispatch(*(*i));
            if (j != j_end) SPEManager::instance()->dispatch(*(*j));
            /// \todo calc scalar parts here
            if (i != i_end) (*i)->wait();
            if (j != j_end) (*j)->wait();
            if (i != i_end) ++i;
            if (j != j_end) ++j;
        }
        for (std::list<SPEInstructionQueue *>::iterator i(iq_lower.begin()), i_end(iq_lower.end()), j(iq_upper.begin()), j_end(iq_upper.end()) ;
                (i != i_end || j != j_end) ; )
        {
            if (i != i_end) delete *i;
            if (j != j_end) delete *j;
            if (i != i_end) ++i;
            if (j != j_end) ++j;
        }
        //bt.take();
        //std::cout<<"wait: "<<bt.sec() - at.sec() << " "<<bt.usec() - at.usec()<<std::endl<<std::endl;

        return result;
    }

    DenseMatrix<float> Product<tags::Cell>::value(const DenseMatrix<float> & a, const DenseMatrix<float> & b)
    {
        CONTEXT("When calculating DenseMatrix<float>-DenseMatrix<float> product (Cell):");

        if (a.columns() != b.rows())
            throw MatrixRowsDoNotMatch(b.rows(), a.columns());

        DenseMatrix<float> result(a.rows(), b.columns());
        std::list<SPEInstruction *> instructions;
        bool use_spe(true);

        union addr
        {
            float * ptr;
            unsigned long long value;
        };

        unsigned a_half_rows = a.rows() / 2;
        for (; ;)
        {
            if ((a_half_rows * b.columns()) % 4 == 0 && (a_half_rows % 4 == 0))
                break;
            else
                a_half_rows++;
        }
        // find a t_size that represents full rows!

        unsigned a_t_size(a.columns() * (4096 / (a.columns())));
        
        while (a_t_size % 4 != 0)
        {
            a_t_size -= a.columns();
        }
        
        unsigned a_div = a_t_size / a.columns(); // number of rows in one transfer
/*
        while(a_div % 4 != 0)
            a_div--;
*/
        //std::cout << "ADIV = " << a_div << std::endl;
        //a_t_size *= 4;
        a_t_size = 4 * (a_div * a.columns());
        //std::cout << a_t_size << std::endl;

        unsigned a_half_size = a_half_rows * a.columns() * 4;
        unsigned a_2nd_half_rows = a.rows() - a_half_rows;
        unsigned a_2nd_half_size = a_2nd_half_rows * a.columns() * 4;
        unsigned a_typed_offset = (a_half_rows * a.columns()) % 4;

        unsigned b_half_rows = b.rows() / 2;
        unsigned b_rest_rows = b.rows() % 2;

        unsigned b_half_cols = ((b.columns() / 2 ) & ~0x3) - 4;
        unsigned ppu_if1_cols = 4;
        unsigned ppu_if2_cols = 4 + (b.columns() & 0x3);
        unsigned b_2nd_half_cols = b.columns() - 4 - b_half_cols - ppu_if2_cols;
/*
        std::cout << "a_half_rows : " << a_half_rows << std::endl;
        std::cout << "a_2nd_half_rows : " << a_2nd_half_rows << std::endl;

        std::cout << "b_half_cols : " << b_half_cols << std::endl;
        std::cout << "b_2nd_half_cols : " << b_2nd_half_cols << std::endl;
        std::cout << "ppu_if1_cols : " << ppu_if1_cols << std::endl;
        std::cout << "ppu_if2_cols : " << ppu_if2_cols << std::endl;
        std::cout << "R0 size: " << a_half_rows * b_half_cols * 4 << std::endl;
*/
        unsigned ppu_columns(ppu_if1_cols + ppu_if2_cols);

        if (a.rows() < 8 || b.columns() < 8)
        {
            use_spe = false;
            ppu_columns = b.columns();
            ppu_if1_cols = b.columns();
            ppu_if2_cols = 0;
            b_half_cols = 0;
        }

        std::vector<SPETransferList> b02_lists = std::vector<SPETransferList>();
        b02_lists.push_back(SPETransferList(2048, 16384));

        std::vector<SPETransferList> b13_lists = std::vector<SPETransferList>();
        b13_lists.push_back(SPETransferList(2048, 16384));

        unsigned b02_nr_lists(0), b13_nr_lists(0);

        if (use_spe)
        {
        for (unsigned i(0) ; i < b.rows() ; ++i)
        {
            addr address1;
            address1.ptr = b.elements() + i * b.columns();
            unsigned off1(address1.value & 0xF);
            address1.value &= ~0xF; // truncate address
            unsigned t_size1 = (off1 == 0 ? b_half_cols * 4 : (b_half_cols * 4) + 16);
            ListElement * retval1(0);
            if ((b02_lists.at(b02_nr_lists).transfer_size() + t_size1) <= 16384)
            {
                retval1 = b02_lists.at(b02_nr_lists).add(address1.ptr, t_size1);
            }

            if (retval1 == 0)
            {
                b02_lists.push_back(SPETransferList(2048, 16384));
                b02_nr_lists++;
                b02_lists.at(b02_nr_lists).add(address1.ptr, t_size1);

            }

            addr address2;
            address2.ptr = ((b.elements() + b_half_cols + ppu_if1_cols) + i * b.columns());
            unsigned off2(address2.value & 0xF);
            address2.value &= ~0xF; // truncate address
            unsigned t_size2 = (off2 == 0 ? b_2nd_half_cols * 4 : (b_2nd_half_cols * 4) + 16);

            ListElement * retval2(0);
            if ((b13_lists.at(b13_nr_lists).transfer_size() + t_size2) <= 16384)
            {
                retval2 = b13_lists.at(b13_nr_lists).add(address2.ptr, t_size2);
            }

            if (retval2 == 0)
            {
                b13_lists.push_back(SPETransferList(2048, 16384));
                b13_nr_lists++;
                b13_lists.at(b13_nr_lists).add(address2.ptr, t_size2);
            }

        }
        }
        unsigned long long b02_sizes[b02_lists.size()] __attribute__((aligned(16)));
        void * b02_ptrs[b02_lists.size()] __attribute__((aligned(16)));
        unsigned long long b02_eahs[b02_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < b02_lists.size() ; i++)
        {
            b02_sizes[i] = b02_lists.at(i).size();
            b02_eahs[i] = b02_lists.at(i).effective_address();
            b02_ptrs[i] = b02_lists.at(i).elements();
        }

        std::vector<SPETransferList> r0_lists = std::vector<SPETransferList>();
        r0_lists.push_back(SPETransferList(2048, 16384));
        std::vector<SPETransferList> r1_lists = std::vector<SPETransferList>();
        r1_lists.push_back(SPETransferList(2048, 16384));
        std::vector<SPETransferList> r2_lists = std::vector<SPETransferList>();
        r2_lists.push_back(SPETransferList(2048, 16384));
        std::vector<SPETransferList> r3_lists = std::vector<SPETransferList>();
        r3_lists.push_back(SPETransferList(2048, 16384));

        unsigned r0_nr_lists(0), r1_nr_lists(0), r2_nr_lists(0), r3_nr_lists(0);
        if (use_spe)
        {
        for (unsigned i(0) ; i < result.rows() ; ++i)
        {
            if (i < a_half_rows)
            {
                addr address1;
                address1.ptr = result.elements() + i * result.columns();
                unsigned off1(address1.value & 0xF);
                address1.value &= ~0xF; // truncate address
                unsigned t_size1 = ((b_half_cols * 4) + 16);

                ListElement * retval1(0);
                if ((r0_lists.at(r0_nr_lists).transfer_size() + t_size1) <= 16384)
                {
                    if (! ((i != 0) && (i % a_div == 0)))
                    {
                        retval1 = r0_lists.at(r0_nr_lists).add(address1.ptr, t_size1);
                    }
                }

                if (retval1 == 0)
                {
                    r0_lists.push_back(SPETransferList(2048, 16384));
                    r0_nr_lists++;
                    r0_lists.at(r0_nr_lists).add(address1.ptr, t_size1);
                }

                addr address2;
                address2.ptr = ((result.elements() + b_half_cols + ppu_if1_cols) + i * result.columns());
                unsigned off2(address2.value & 0xF);
                address2.value &= ~0xF; // truncate address
                unsigned t_size2 = (b_2nd_half_cols * 4) + 16;

                ListElement * retval2(0);
                if ((r1_lists.at(r1_nr_lists).transfer_size() + t_size2) <= 16384)
                {
                    if (! ((i != 0) && (i % a_div == 0)))
                    {
                        retval2 = r1_lists.at(r1_nr_lists).add(address2.ptr, t_size2);
                    }
                }

                if (retval2 == 0)
                {
                    r1_lists.push_back(SPETransferList(2048, 16384));
                    r1_nr_lists++;
                    r1_lists.at(r1_nr_lists).add(address2.ptr, t_size2);
                }
            }
            else
            {
                addr address1;
                address1.ptr = result.elements() + i * result.columns();
                unsigned off1(address1.value & 0xF);
                address1.value &= ~0xF; // truncate address
                unsigned t_size1 = (b_half_cols * 4) + 16;
                ListElement * retval1(0);
                if ((r2_lists.at(r2_nr_lists).transfer_size() + t_size1) <= 16384)
                {
                    if (! (((i - a_half_rows) != 0) && ((i - a_half_rows) % a_div == 0)))
                    {
                        retval1 = r2_lists.at(r2_nr_lists).add(address1.ptr, t_size1);
                    }
                }

                if (retval1 == 0)
                {
                    r2_lists.push_back(SPETransferList(2048, 16384));
                    r2_nr_lists++;
                    r2_lists.at(r2_nr_lists).add(address1.ptr, t_size1);
                }

                addr address2;
                address2.ptr = ((result.elements() + b_half_cols + ppu_if1_cols) + i * result.columns());
                unsigned off2(address2.value & 0xF);
                address2.value &= ~0xF; // truncate address
//                unsigned t_size2 = (off2 == 0 ? b_2nd_half_cols * 4 : (b_2nd_half_cols * 4) + 16);
                unsigned t_size2 = (b_2nd_half_cols * 4) + 16;

                ListElement * retval2(0);
                if ((r3_lists.at(r3_nr_lists).transfer_size() + t_size2) <= 16384)
                {
                    if (! (((i - a_half_rows) != 0) && ((i - a_half_rows) % a_div == 0)))
                    {
                        retval2 = r3_lists.at(r3_nr_lists).add(address2.ptr, t_size2);
                    }
                }

                if (retval2 == 0)
                {
                    r3_lists.push_back(SPETransferList(2048, 16384));
                    r3_nr_lists++;
                    r3_lists.at(r3_nr_lists).add(address2.ptr, t_size2);
                }
            }

        }
        }
        unsigned long long r0_sizes[r0_lists.size()] __attribute__((aligned(16)));
        void * r0_ptrs[r0_lists.size()] __attribute__((aligned(16)));
        unsigned long long r0_eahs[r0_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < r0_lists.size() ; i++)
        {
            r0_sizes[i] = r0_lists.at(i).size();
            r0_eahs[i] = r0_lists.at(i).effective_address();
            r0_ptrs[i] = r0_lists.at(i).elements();
        }

        unsigned long long r1_sizes[r1_lists.size()] __attribute__((aligned(16)));
        unsigned long long r1_eahs[r1_lists.size()] __attribute__((aligned(16)));
        void * r1_ptrs[r1_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < r1_lists.size() ; i++)
        {
            r1_sizes[i] = r1_lists.at(i).size();
            r1_eahs[i] = r1_lists.at(i).effective_address();
            r1_ptrs[i] = r1_lists.at(i).elements();
        }
        unsigned long long r2_sizes[r2_lists.size()] __attribute__((aligned(16)));
        void * r2_ptrs[r2_lists.size()] __attribute__((aligned(16)));
        unsigned long long r2_eahs[r2_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < r2_lists.size() ; i++)
        {
            r2_sizes[i] = r2_lists.at(i).size();
            r2_eahs[i] = r2_lists.at(i).effective_address();
            r2_ptrs[i] = r2_lists.at(i).elements();
        }

        unsigned long long r3_sizes[r3_lists.size()] __attribute__((aligned(16)));
        unsigned long long r3_eahs[r3_lists.size()] __attribute__((aligned(16)));
        void * r3_ptrs[r3_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < r3_lists.size() ; i++)
        {
            r3_sizes[i] = r3_lists.at(i).size();
            r3_eahs[i] = r3_lists.at(i).effective_address();
            r3_ptrs[i] = r3_lists.at(i).elements();
        }

        // Shared Operands
        Operand oa0 = { a.elements() };

        Operand ob0 = { a_half_size / a_t_size }; // # of DBs
        Operand oc0 = { a_half_size % a_t_size }; // rest

        if (oc0.u == 0)
        {
            if (ob0.u > 0)
                oc0.u = a_t_size;
        }
        else
        {
            ob0.u++;
        }
        //std::cout << "TRANSFERS FOR A: " << ob0.u << std::endl;
        Operand oa1 = { a.elements() + (a_half_rows * a.columns()) }; //truncated
        oa1.u &= ~0xF;

        Operand ob1 = { (a_2nd_half_size + (4 * a_typed_offset)) / a_t_size  }; // # of DBs
        Operand oc1 = { (a_2nd_half_size + (4 * a_typed_offset)) % a_t_size  }; // rest (m_of_s on SPU)

        if (oc1.u == 0)
        {
            if (ob1.u > 0)
                oc1.u = a_t_size;
        }
        else
        {
            ob1.u++;
        }

        Operand oi = { a.columns() };

        { // DISPATCH FOR R[0] and R[2]
            Operand od = { &b02_ptrs };
            Operand oe = { &b02_sizes };
            Operand of = { &b02_eahs };
            Operand og = { b02_lists.size() };
            //std::cout << "B02 LIST SIZE: " << og.u << std::endl;
            Operand oh = { b.elements() + b.columns() };
            oh.u &= 0xF; // Want the offset of the first row!

            Operand oj = { b_half_cols };

            Operand ok0 = { r0_lists.size() };
            //std::cout << "R0 LIST SIZE: " << ok0.u << std::endl;

            Operand ok1 = { r2_lists.size() };

            Operand ol0 = { &r0_ptrs };
            Operand ol1 = { &r2_ptrs };

            Operand om0 = { &r0_sizes };
            Operand om1 = { &r2_sizes };

            Operand on0 = { &r0_eahs };
            Operand on1 = { &r2_eahs };

            SPEInstruction * instruction0 = new SPEInstruction(oc_product_dense_matrix_dense_matrix_float, a_t_size, oa0, ob0, oc0, od, oe, of, og, oh, oi,
                    oj, ok0, ol0, om0, on0);
            SPEInstruction * instruction1 = new SPEInstruction(oc_product_dense_matrix_dense_matrix_float, a_t_size, oa1, ob1, oc1, od, oe, of, og, oh, oi,
                    oj, ok1, ol1, om1, on1);

            if (use_spe)
            {
            SPEManager::instance()->dispatch(*instruction0);
            instructions.push_back(instruction0);
            SPEManager::instance()->dispatch(*instruction1);
            instructions.push_back(instruction1);
            }
        }

        unsigned long long b13_sizes[b13_lists.size()] __attribute__((aligned(16)));
        unsigned long long b13_eahs[b13_lists.size()] __attribute__((aligned(16)));
        void * b13_ptrs[b13_lists.size()] __attribute__((aligned(16)));

        for (unsigned i(0) ; i < b13_lists.size() ; i++)
        {
            b13_sizes[i] = b13_lists.at(i).size();
            b13_eahs[i] = b13_lists.at(i).effective_address();
            b13_ptrs[i] = b13_lists.at(i).elements();
        }

        { // DISPATCH FOR R[1] and R[3]

            Operand od = { &b13_ptrs };
            Operand oe = { &b13_sizes };
            Operand of = { &b13_eahs };

            Operand og = { b13_lists.size() };

            Operand oh = { (b.elements() + b_half_cols + ppu_if1_cols) + b.columns() };
            oh.u &= 0xF;

            Operand oj = { b_2nd_half_cols };

            Operand ok0 = { r1_lists.size() };
            Operand ok1 = { r3_lists.size() };

            Operand ol0 = { &r1_ptrs };
            Operand ol1 = { &r3_ptrs };

            Operand om0 = { &r1_sizes };
            Operand om1 = { &r3_sizes };

            Operand on0 = { &r1_eahs };
            Operand on1 = { &r3_eahs };

            SPEInstruction * instruction2 = new SPEInstruction(oc_product_dense_matrix_dense_matrix_float, a_t_size, oa0, ob0, oc0, od, oe, of, og, oh, oi, oj,
                    ok0, ol0, om0, on0);
            SPEInstruction * instruction3 = new SPEInstruction(oc_product_dense_matrix_dense_matrix_float, a_t_size, oa1, ob1, oc1, od, oe, of, og, oh, oi, oj,
                    ok1, ol1, om1, on1);

            if (use_spe)
            {
            SPEManager::instance()->dispatch(*instruction2);
            instructions.push_back(instruction2);
            SPEManager::instance()->dispatch(*instruction3);
            instructions.push_back(instruction3);
            }
        }

        float * b_if1_start;
        float * b_if2_start;

        float cols[ppu_columns][result.rows()];

        for (unsigned x(0) ; x < ppu_columns ; ++x)
        {
            TypeTraits<float>::fill(cols[x], result.rows(), 0.0f);
        }

        float * a_elem = a.elements();
        //for(Matrix<float>::ConstElementIterator j(a.begin_elements()), j_end(a.end_elements()) ; j != j_end ; ++j)
        for(unsigned j(0) ; j < a.rows() * a.columns() ; ++j)
        {
            if (j % a.columns() == 0)
            {
                // Reset in case of new Row in A
                b_if1_start = b.elements() + b_half_cols;
                b_if2_start = b.elements() + (b.columns() - ppu_if2_cols);
            }

            unsigned i(0);
            for( ; i < ppu_if1_cols ; i++)
            {
                cols[i][j / a.columns()] += *a_elem * b_if1_start[i];
            }

            for(unsigned k(0) ; k < ppu_if2_cols ; k++)
            {
                cols[i+k][j / a.columns()] += *a_elem * b_if2_start[k];
            }

            b_if1_start += b.columns();
            b_if2_start += b.columns();
            a_elem++;

        }

        for(std::list<SPEInstruction *>::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; ++i)
        {
                (*i)->wait();

                delete *i;
        }

        // COPY cols to Matrix.
        unsigned i(0);
        for (unsigned j(b_half_cols) ; i < ppu_if1_cols ; i++, j++)
        {
            DenseVectorSlice<float> slice(result.column(j));
            for(DenseVectorSlice<float>::ElementIterator ei(slice.begin_elements()), ei_end(slice.end_elements()) ; ei != ei_end ; ++ei)
            {
                *ei = cols[i][ei.index()];
            }
        }
        for (unsigned k(0), j(b.columns() - ppu_if2_cols) ; k < ppu_if2_cols ; k++, j++)
        {
            DenseVectorSlice<float> slice(result.column(j));
            for(DenseVectorSlice<float>::ElementIterator ei(slice.begin_elements()), ei_end(slice.end_elements()) ; ei != ei_end ; ++ei)
            {
                *ei = cols[i+k][ei.index()];
            }
        }

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
