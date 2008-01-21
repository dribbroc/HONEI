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
        CONTEXT("When calculating node distances for DenseMatrix<float> (Cell):");

        DenseMatrix<float> result(pos_matrix.rows(), pos_matrix.rows(), float(0));

        for(unsigned i(0) ; i < pos_matrix.rows() ; i++)
        {
            unsigned off(pos_matrix.rows() % 16);
            DenseVector<float> result_row(pos_matrix.rows() + 16 - off, float(0));

            Operand oa = { pos_matrix.elements() };

            Operand ob, oc;
            ob.u = (pos_matrix.rows() * 2) / (16384 / sizeof(float));
            oc.u = (pos_matrix.rows() * 2) % (16384 / sizeof(float));
            oc.u *= 4;

            Operand od, oe;
            od.f = *(pos_matrix.elements() + (i * 2));
            oe.f = *(pos_matrix.elements() + (i * 2) + 1);

            Operand of = { result_row.elements() };
            Operand og = { result_row.size() };

            bool use_spe(true);
            if (0 == oc.u)
            {
                if (ob.u > 0)
                {
                    oc.u = 16 * 1024;
                }
                else
                {
                    use_spe = false;
                }
            }
            else
            {
                ++ob.u;
            }

            if (use_spe)
            {
                SPEInstruction instruction(oc_node_distance_float, 16384, oa, ob, oc, od, oe, of, og);

                SPEManager::instance()->dispatch(instruction);

                instruction.wait();
            }

            float * address = result.elements() + (result.rows() * i);
            DenseVectorRange<float> res(result_row, pos_matrix.rows(), 0);
            Vector<float>::ElementIterator i(res.begin_elements());

            for (int j(0); j < result.rows() ; j++)
            {
                address[j] = *i;
                ++i;
            }
        }

        return result;
    }

    void
    NodeDistance<tags::Cell>::value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<bool> & neighbours,
        DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
    {
        CONTEXT("When calculating node distances for DenseMatrix<float> (Cell):");

        MutableMatrix<float>::ElementIterator e(inv_square_dist.begin_elements());
        MutableMatrix<float>::ElementIterator f(square_dist.begin_elements());
        Matrix<bool>::ConstElementIterator g(neighbours.begin_elements());
        float square_force_range(repulsive_force_range * repulsive_force_range);

        DenseMatrix<float> result(pos_matrix.rows(), pos_matrix.rows(), float(0));

        for(unsigned i(0) ; i < pos_matrix.rows() ; i++)
        {
            unsigned off(pos_matrix.rows() % 16);
            DenseVector<float> result_row(pos_matrix.rows() + 16 - off, float(0));

            Operand oa = { pos_matrix.elements() };

            Operand ob, oc;
            ob.u = (pos_matrix.rows() * 2) / (16384 / sizeof(float));
            oc.u = (pos_matrix.rows() * 2) % (16384 / sizeof(float));
            oc.u *= 4;

            Operand od, oe;
            od.f = *(pos_matrix.elements() + (i * 2));
            oe.f = *(pos_matrix.elements() + (i * 2) + 1);

            Operand of = { result_row.elements() };
            Operand og = { result_row.size() };

            bool use_spe(true);
            if (0 == oc.u)
            {
                if (ob.u > 0)
                {
                    oc.u = 16 * 1024;
                }
                else
                {
                    use_spe = false;
                }
            }
            else
            {
                ++ob.u;
            }

            if (use_spe)
            {
                SPEInstruction instruction(oc_node_distance_float, 16384, oa, ob, oc, od, oe, of, og);

                SPEManager::instance()->dispatch(instruction);

                instruction.wait();
            }

            float * address = result.elements() + (result.rows() * i);
            DenseVectorRange<float> res(result_row, pos_matrix.rows(), 0);
            Vector<float>::ElementIterator i(res.begin_elements());

            for (int j(0); j < result.rows() ; j++)
            {
                address[j] = *i;
                if (*i < square_force_range && *i > std::numeric_limits<float>::epsilon())
                {
                    *e = 1 / *i;
                }
                else
                {
                    *e = 0.0f;
                }

                if (*g)
                {
                    *f = *i;
                }

                ++i; ++e; ++f; ++g;
            }

        }
    }

    void
    NodeDistance<tags::Cell>::value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<float> & edge_weights,
        DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
    {
        CONTEXT("When calculating node distances for DenseMatrix<float> (Cell):");

        MutableMatrix<float>::ElementIterator e(inv_square_dist.begin_elements());
        MutableMatrix<float>::ElementIterator f(square_dist.begin_elements());
        Matrix<float>::ConstElementIterator g(edge_weights.begin_elements());
        float square_force_range(repulsive_force_range * repulsive_force_range);

        DenseMatrix<float> result(pos_matrix.rows(), pos_matrix.rows(), float(0));

        for(unsigned i(0) ; i < pos_matrix.rows() ; i++)
        {
            unsigned off(pos_matrix.rows() % 16);
            DenseVector<float> result_row(pos_matrix.rows() + 16 - off, float(0));

            Operand oa = { pos_matrix.elements() };

            Operand ob, oc;
            ob.u = (pos_matrix.rows() * 2) / (16384 / sizeof(float));
            oc.u = (pos_matrix.rows() * 2) % (16384 / sizeof(float));
            oc.u *= 4;

            Operand od, oe;
            od.f = *(pos_matrix.elements() + (i * 2));
            oe.f = *(pos_matrix.elements() + (i * 2) + 1);

            Operand of = { result_row.elements() };
            Operand og = { result_row.size() };

            bool use_spe(true);
            if (0 == oc.u)
            {
                if (ob.u > 0)
                {
                    oc.u = 16 * 1024;
                }
                else
                {
                    use_spe = false;
                }
            }
            else
            {
                ++ob.u;
            }

            if (use_spe)
            {
                SPEInstruction instruction(oc_node_distance_float, 16384, oa, ob, oc, od, oe, of, og);

                SPEManager::instance()->dispatch(instruction);

                instruction.wait();
            }

            float * address = result.elements() + (result.rows() * i);
            DenseVectorRange<float> res(result_row, pos_matrix.rows(), 0);
            Vector<float>::ElementIterator i(res.begin_elements());

            for (int j(0); j < result.rows() ; j++)
            {
                address[j] = *i;

                if (*i < square_force_range && *i > std::numeric_limits<float>::epsilon())
                {
                    *e = 1 / *i;
                }
                else
                {
                    *e = 0.0f;
                }

                if (*g > std::numeric_limits<float>::epsilon())
                {
                    *f = *i;
                }

                ++i; ++e; ++f; ++g;
            }

        }
    }


}
