/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/libgraph/node_distance.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_instruction.hh>
#include <honei/libutil/spe_manager.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float>
    NodeDistance<tags::Cell>::value(const DenseMatrix<float> & pos_matrix)
    {
        CONTEXT("When calculating node distances (WKK) for DenseMatrix<float> (Cell):");

        std::list<SPEInstruction *> instructions;
        DenseMatrix<float> result(pos_matrix.rows(), pos_matrix.rows());

        unsigned a_rows = 16384 / (pos_matrix.rows() * 4);
        if (a_rows > pos_matrix.rows())
            a_rows = pos_matrix.rows();

        while ((a_rows % 2) != 0 || (a_rows * pos_matrix.rows() % 4 != 0))
            a_rows--;

        unsigned needed_iterations = pos_matrix.rows() / a_rows;
        unsigned ppu_start_row = needed_iterations * a_rows;
        unsigned ppu_rows = (pos_matrix.rows() - ppu_start_row);

        for (unsigned i(0) ; i < needed_iterations ; i++)
        {
            Operand oa = { pos_matrix.elements() + (i * a_rows * 2) };
            Operand ob = { pos_matrix.elements() };
            Operand oc;
            oc.u = a_rows * 2;

            Operand od, oe; // Transfer sizes for B = pos_matrix (full matrix with partial transfers)
            od.u = ((pos_matrix.rows() * 2)) / (16384 / sizeof(float));
            oe.u = ((pos_matrix.rows() * 2)) % (16384 / sizeof(float));

            oc.u *= 4;
            oe.u *= 4;
            Operand of = { result.elements() + (i * a_rows * pos_matrix.rows()) };
            Operand og;
            og.u = pos_matrix.rows();

            bool use_spe(true);

            if (0 == oc.u)
            {
                if (pos_matrix.rows() > 0)
                {
                    oc.u = 16 * 1024;
                }
                else
                {
                    use_spe = false;
                }
            }

            if (0 == oe.u)
            {
                if (od.u > 0)
                {
                    oe.u = 16 * 1024;
                }
            }
            else
            {
                ++od.u;
            }

            if (use_spe && i % 6 != 5)
            {
                SPEInstruction * instruction = new SPEInstruction(oc_node_distance_float_wkk, 16384, oa, ob, oc, od, oe, of, og);

                SPEManager::instance()->dispatch(*instruction);

                instructions.push_back(instruction);
            }
            else
            {
                for(std::list<SPEInstruction * >::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; ++i)
                {
                        (*i)->wait();
                }

                SPEInstruction * instruction = new SPEInstruction(oc_node_distance_float_wkk, 16384, oa, ob, oc, od, oe, of, og);

                SPEManager::instance()->dispatch(*instruction);

                instructions.push_back(instruction);

            }

        }

        unsigned starting_point(ppu_start_row * pos_matrix.rows());

        MutableMatrix<float>::ElementIterator j(result.element_at(starting_point));
        for(Matrix<float>::ConstElementIterator k(pos_matrix.element_at(ppu_start_row * 2)), k_end(pos_matrix.end_elements()) ; k != k_end ; )
        {
            float a_x = *k;
            ++k;
            float a_y = *k;
            ++k;

            for (Matrix<float>::ConstElementIterator i(pos_matrix.begin_elements()), i_end(pos_matrix.end_elements()) ; i != i_end ; )
            {
                float x_temp = a_x - *i;
                ++i;
                float y_temp = a_y - *i;
                ++i;
                x_temp *= x_temp;
                y_temp *= y_temp;

                *j = x_temp + y_temp;

                ++j;
            }
        }

        for(std::list<SPEInstruction * >::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; ++i)
        {
                (*i)->wait();
                delete *i;
        }

        return result;
    }

    void
    NodeDistance<tags::Cell>::value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<float> & edge_weights,
        DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
    {
        CONTEXT("When calculating node distances (WFR) for DenseMatrix<float> (Cell):");

        std::list<SPEInstruction *> instructions;
        float square_force_range(repulsive_force_range * repulsive_force_range);

        unsigned a_rows = 16384 / (pos_matrix.rows() * 4);
        if (a_rows > pos_matrix.rows())
            a_rows = pos_matrix.rows();

        while ((a_rows % 2) != 0 || (a_rows * pos_matrix.rows() % 4 != 0))
            a_rows--;

        unsigned needed_iterations = pos_matrix.rows() / a_rows;
        unsigned ppu_start_row = needed_iterations * a_rows;
        unsigned ppu_rows = (pos_matrix.rows() - ppu_start_row);

        for (unsigned i(0) ; i < needed_iterations ; i++)
        {
            Operand oa = { pos_matrix.elements() + (i * a_rows * 2) };
            Operand ob = { pos_matrix.elements() };
            Operand oc;
            oc.u = a_rows * 2;

            Operand od, oe; // Transfer sizes for B = pos_matrix (full matrix with partial transfers)
            od.u = ((pos_matrix.rows() * 2)) / (16384 / sizeof(float));
            oe.u = ((pos_matrix.rows() * 2)) % (16384 / sizeof(float));

            oc.u *= 4;
            oe.u *= 4;

            Operand of = { inv_square_dist.elements() + (i * a_rows * pos_matrix.rows()) };
            Operand og = { square_dist.elements() + (i * a_rows * pos_matrix.rows()) };
            Operand oh = { edge_weights.elements() + (i * a_rows * pos_matrix.rows()) };
            Operand oi, oj, ok;
            oi.u = pos_matrix.rows();
            oj.f = std::numeric_limits<float>::epsilon();
            ok.f = square_force_range;

            bool use_spe(true);

            if (0 == oc.u)
            {
                if (pos_matrix.rows() > 0)
                {
                    oc.u = 16 * 1024;
                }
                else
                {
                    use_spe = false;
                }
            }

            if (0 == oe.u)
            {
                if (od.u > 0)
                {
                    oe.u = 16 * 1024;
                }
            }
            else
            {
                ++od.u;
            }

            if (use_spe && i % 6 != 5)
            {
                SPEInstruction * instruction = new SPEInstruction(oc_node_distance_float_wfr, 16384, oa, ob, oc, od, oe, of, og, oh, oi, oj, ok);

                SPEManager::instance()->dispatch(*instruction);

                instructions.push_back(instruction);
            }
            else
            {
                for(std::list<SPEInstruction * >::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; ++i)
                {
                        (*i)->wait();
                }

                SPEInstruction * instruction = new SPEInstruction(oc_node_distance_float_wfr, 16384, oa, ob, oc, od, oe, of, og, oh, oi, oj, ok);

                SPEManager::instance()->dispatch(*instruction);

                instructions.push_back(instruction);
            }
        }

        DenseMatrix<float> tiny_result(ppu_rows, pos_matrix.rows());

        unsigned starting_point(ppu_start_row * pos_matrix.rows());
        MutableMatrix<float>::ElementIterator e(inv_square_dist.element_at(starting_point));
        MutableMatrix<float>::ElementIterator f(square_dist.element_at(starting_point));
        Matrix<float>::ConstElementIterator g(edge_weights.element_at(starting_point));

        MutableMatrix<float>::ElementIterator j(tiny_result.begin_elements());
        for(Matrix<float>::ConstElementIterator k(pos_matrix.element_at(ppu_start_row * 2)), k_end(pos_matrix.end_elements()) ; k != k_end ; )
        {
            float a_x = *k;
            ++k;
            float a_y = *k;
            ++k;

            for (Matrix<float>::ConstElementIterator i(pos_matrix.begin_elements()), i_end(pos_matrix.end_elements()) ; i != i_end ; )
            {
                float x_temp = a_x - *i;
                ++i;
                float y_temp = a_y - *i;
                ++i;
                x_temp *= x_temp;
                y_temp *= y_temp;

                *j = x_temp + y_temp;

                if (*j < square_force_range && *j > std::numeric_limits<float>::epsilon())
                {
                    *e = 1 / *j;
                }
                else
                {
                    *e = 0.0f;
                }

                if (*g > std::numeric_limits<float>::epsilon())
                {
                    *f = *j;
                }

                ++j; ++e; ++f; ++g;
            }
        }

        for(std::list<SPEInstruction * >::iterator i(instructions.begin()), i_end(instructions.end()) ; i != i_end ; ++i)
        {
                (*i)->wait();
                delete *i;
        }
    }
}
