/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <honei/lbm/solver_labswe_grid.hh>

using namespace honei::lbm;
using namespace honei::lbm::lbm_boundary_types;
namespace honei
{

    template<typename ResPrec_>
    SolverLABSWEGrid<tags::CPU::MultiCore, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>::SolverLABSWEGrid(PackedGridData<D2Q9, ResPrec_> * data, PackedGridInfo<D2Q9> * info, ResPrec_ dx, ResPrec_ dy, ResPrec_ dt) :
        _parts(4), /// \todo use Configuration
        _data(data),
        _info(info)
    {
        CONTEXT("When creating LABSWE solver:");
        GridPartitioner<D2Q9, ResPrec_>::decompose(_parts, *_info, *_data, _info_list, _data_list, _fringe_list);

        for(unsigned long i(0) ; i < _parts ; ++i)
        {
            _solver_list.push_back(new SolverLABSWEGrid<tags::CPU, ResPrec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>(&_data_list[i], &_info_list[i], 1., 1., 1.));
        }
    }

    template<typename ResPrec_>
    SolverLABSWEGrid<tags::CPU::MultiCore, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>::~SolverLABSWEGrid()
        {
            CONTEXT("When destroying LABSWE solver.");
        }

    template<typename ResPrec_>
    void
    SolverLABSWEGrid<tags::CPU::MultiCore, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>::do_preprocessing()
    {
        CONTEXT("When performing LABSWE preprocessing.");
        for (unsigned long i(0) ; i < _parts ; ++i)
        {
            _solver_list.at(i)->do_preprocessing();
        }
        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
    }

    template<typename ResPrec_>
    void
    SolverLABSWEGrid<tags::CPU::MultiCore, ResPrec_, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>::solve()
    {
        for (unsigned long i(0) ; i < _parts ; ++i)
        {
            _solver_list.at(i)->solve();
        }
        GridPartitioner<D2Q9, ResPrec_>::synch(*_info, *_data, _info_list, _data_list, _fringe_list);
        /// \todo remove compose - it is only necessary if one must read the data
        GridPartitioner<D2Q9, ResPrec_>::compose(*_info, *_data, _info_list, _data_list);
    }

    template
    class SolverLABSWEGrid<tags::CPU::MultiCore, float, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>;

    template
    class SolverLABSWEGrid<tags::CPU::MultiCore, double, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>;
}
