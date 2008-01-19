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

#include <cell/cell.hh>
#include <libgraph/node_distance.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace honei
{
    using namespace cell;

    DenseMatrix<float>
    NodeDistance<tags::Cell>::value(const DenseMatrix<float> & a)
    {
        CONTEXT("When calculating node distances for DenseMatrix<float> (Cell):");

        DenseMatrix<float> result(a.rows(), a.rows(), float(0));

        for(unsigned i(0) ; i < a.rows() ; i++)
        {
            unsigned off(a.rows() % 16);
            DenseVector<float> result_row(a.rows() + 16 - off, float(0));

            Operand oa = { a.elements() };

            Operand ob, oc;
            ob.u = (a.rows() * 2) / (16384 / sizeof(float));
            oc.u = (a.rows() * 2) % (16384 / sizeof(float));
            oc.u *= 4;

            Operand od, oe;
            od.f = *(a.elements() + (i * 2));
            oe.f = *(a.elements() + (i * 2) + 1);

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
            DenseVectorRange<float> res(result_row, a.rows(), 0);
            Vector<float>::ElementIterator i(res.begin_elements());

            for (int j(0); j < result.rows() ; j++)
            {
                std::cout << *i << std::endl;
                address[j] = *i;
                ++i;
            }
        }

        return result;
    }
}
