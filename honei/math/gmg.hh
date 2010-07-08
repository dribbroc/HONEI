/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.tu-dortmund.de>
 *
 * This file is part of the HONEI math C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MATH_GUARD_GMG_HH
#define MATH_GUARD_GMG_HH 1

#include <vector>
#include <tr1/functional>
#include <honei/la/dense_vector.hh>
#include <honei/math/methods.hh>


namespace honei
{
    template<typename Prec_, typename MatrixType_, typename ProlMatrixType_>
    struct GMGInfo
    {
        public:
            std::vector<unsigned long> cycle;
            std::vector<std::tr1::function<void ()> > smoother_functors;
            std::vector<std::tr1::function<void ()> > transfer_functors;
            std::vector<std::tr1::function<void ()> > coarse_solver_functors;

            std::vector<MatrixType_> system_matrices;
            std::vector<ProlMatrixType_> prolongation_matrices;
            std::vector<ProlMatrixType_> restriction_matrices;
            std::vector<DenseVector<Prec_> > rhs_vectors;
            std::vector<DenseVector<Prec_> > result_vectors;
            std::vector<DenseVector<Prec_> > c_vectors;
            std::vector<DenseVector<Prec_> > d_vectors;

            std::vector<MatrixType_> precon_matrices;
            std::vector<DenseVector<Prec_> > precon_vectors;

            unsigned long start_level;
            unsigned long min_level;
            unsigned long max_level;
            unsigned long end_level;
            unsigned long max_iters_global;
            unsigned long max_iters_smoother;
            unsigned long max_iters_coarse_solver;
            Prec_ tolerance;
            Prec_ tolerance_coarse;
            Prec_ adaptive_correction_factor;
    };

    template<typename Prec_,
             typename MatrixType_,
             typename ProlMatrixType_,
             typename SmootherType_,
             typename PreconType_,
             typename CoarseSolverType_,
             typename ProlType_,
             typename ResType_,
             typename NormType_>
    struct GMGInfoFactory
    {
        public:
            static GMGInfo<Prec_,
                           MatrixType_,
                           ProlMatrixType_> create(std::tr1::function<void (std::vector<MatrixType_> &,
                                                                            std::vector<ProlMatrixType_> &,
                                                                            std::vector<ProlMatrixType_> &,
                                                                            std::vector<DenseVector<Prec_> > &,
                                                                            unsigned long, unsigned long)> problem_factory,
                                                   unsigned long min_level,
                                                   unsigned long max_level,
                                                   unsigned long start_level,
                                                   unsigned long end_level,
                                                   unsigned long max_iters_global,
                                                   unsigned long max_iters_smoother,
                                                   unsigned long max_iters_coarse_solver,
                                                   Prec_ tolerance,
                                                   Prec_ tolerance_coarse,
                                                   Prec_ adaptive_correction_factor,
                                                   std::vector<unsigned long> & cycle)
            {
                GMGInfo<Prec_, MatrixType_, ProlMatrixType_> info;
                problem_factory(info.system_matrices, info.prolongation_matrices, info.restriction_matrices, info.rhs_vectors, min_level, max_level);
                info.cycle = cycle;

                for (unsigned long i(0) ; i < max_level - min_level ; ++i)
                {
                    info.precon_matrices.push_back(PreconType_::value(info.system_matrices.at(i)));
                }
                return info;
            }
    };

    template<typename Prec_,
             typename MatrixType_,
             typename ProlMatrixType_>
    class GMGState
    {
        private:
            GMGInfo<Prec_, MatrixType_, ProlMatrixType_> & _info;

        public:
            unsigned long current_level;
            GMGState(GMGInfo<Prec_, MatrixType_, ProlMatrixType_> & info) :
                _info(info),
                current_level(info.start_level)
            {
            }

            void descent(std::tr1::function<void ()> & smoother_functor, std::tr1::function<void ()> & transfer_functor)
            {
                smoother_functor();
                transfer_functor();
            }

            void rise(std::tr1::function<void ()> & smoother_functor, std::tr1::function<void ()> & transfer_functor)
            {
                transfer_functor();
                smoother_functor();
            }

            void solve_coarse(std::tr1::function<void ()> & coarse_smoother_functor)
            {
                ///TODO:Coarse Grid corrector?
                coarse_smoother_functor();
            }
    };

    template<typename Tag_, typename Mode_>
    struct GMG
    {
        public:
            template<typename Prec_, typename MatrixType_, typename ProlMatrixType_>
                static void value(GMGInfo<Prec_, MatrixType_, ProlMatrixType_> & info)
                {
                    GMGState<Prec_, MatrixType_, ProlMatrixType_> state(info);
                    std::vector<unsigned long>::iterator cycle_iterator(info.cycle.begin());
                    std::vector<std::tr1::function<void ()> >::iterator smoother_functors_iterator(info.smoother_functors.begin());
                    std::vector<std::tr1::function<void ()> >::iterator coarse_solver_functors_iterator(info.coarse_solver_functors.begin());
                    std::vector<std::tr1::function<void ()> >::iterator transfer_functors_iterator(info.transfer_functors.begin());
                    for(; cycle_iterator != info.cycle.end() ; ++cycle_iterator, ++smoother_functors_iterator, ++transfer_functors_iterator)
                    {
                        if(*cycle_iterator < state.current_level)
                            state.descent(*smoother_functors_iterator, *transfer_functors_iterator);
                        else
                            state.rise(*smoother_functors_iterator, *transfer_functors_iterator);

                        if(state.current_level == info.min_level)
                        {
                            state.solve_coarse(*coarse_solver_functors_iterator);
                            ++coarse_solver_functors_iterator;
                        }

                    }
                }
    };
}

#endif
