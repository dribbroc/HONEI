/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBLBM_GUARD_SOLVER_LABSWE_HH
#define LIBLBM_GUARD_SOLVER_LABSWE_HH 1

/**
 * \file
 * Implementation of a SWE solver using LBM.
 *
 * \ingroup grpliblbm
 **/

#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>

namespace honei
{

	namespace lbm_lattice_types
	{
		class D2Q9;
	}
	namespace lbm_grid_types
	{
		class RECTANGULAR;
	}
	namespace lbm_source_types
	{
		class SIMPLE;
	}
	namespace lbm_source_shemes
	{
		class BASIC;
		class CENTERED;
	}
	namespace lbm_boundary_types
	{
		class NOSLIP_PERIODIC;
	}

	template<typename Tag_,
			 typename ResPrec_,
			 typename SourceType_,
			 typename SourceSheme_,
			 typename GridType_,
			 typename LatticeType_,
			 typename BoundaryType_>
	class SolverLABSWE
	{
	};

	template<typename Tag_, typename ResPrec_>
	class SolverLABSWE<Tag_, ResPrec_, lbm_source_types::SIMPLE, lbm_source_shemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC>
	{
		private:
			/** Global variables.
			  *
			 **/
			ResPrec_ _relaxation_time, _delta_x, _delta_y, _delta_t;
			DenseMatrix<ResPrec_>* _height;
			DenseMatrix<ResPrec_>* _u;
			DenseMatrix<ResPrec_>* _v;

			unsigned long _grid_width, _grid_height;

			DenseVector<ResPrec_>* _equilibrium_distribution;
			DenseVector<ResPrec_>* _source_x;
			DenseVector<ResPrec_>* _source_y;

			/** Global constants.
			  *
			 **/
			ResPrec_ _n_alpha, _e, _gravity;
			DenseVector<ResPrec_>* _distribution_vector_x;
			DenseVector<ResPrec_>* _distribution_vector_y;


		public:
			SolverLABSWE(ResPrec_ dx, ResPrec_ dy, ResPrec_ dt, unsigned long gx, unsigned long gy)
			{

			}
	};

}
#endif
