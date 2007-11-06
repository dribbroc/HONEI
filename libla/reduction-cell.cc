/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Till Barz <till.barz@uni-dortmund.de>
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
#include <libla/reduction.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_manager.hh>

namespace honei
{
    float 
    Reduction<rt_sum,tags::Cell>::value(const DenseVector<float> & a)
    {
        CONTEXT("When reducing DenseVector<float> to Scalar by Sum (Cell):");

        float result;
        
        Operand oa = { a.elements() };
        Operand oc = { &result };
        
        SPEInstruction instruction(oc_dense_float_vector_reduction_sum, a.size(), oa, oc);

        SPEManager::instance()->dispatch(instruction);

        instruction.wait();
       
        return result;
    }
    


}
