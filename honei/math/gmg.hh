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

#include<list>
#include<tr1/functional>
#include<honei/math/methods.hh>


namespace honei
{
    template<typename Prec_>
    struct GMGInfo
    {
        public:
            std::list<unsigned long> cycle;
            std::list<std::tr1::function<void ()> > smoother_functors;
            std::list<std::tr1::function<void ()> > transfer_functors;

            std::list<std::tr1::function<void ()> > coarse_solver_functors;
            unsigned long start_level;
            unsigned long min_level;
    };

    template<typename Prec_, typename CycleType_>
    struct GMGInfoFactory;

    template<typename Prec_>
    struct GMGInfoFactory<Prec_, methods::CYCLE::V>
    {
        static GMGInfo<Prec_> create()
        {
            GMGInfo<Prec_> result;
            return result;
        }
    };

    template<typename Prec_>
    class GMGState
    {
        private:
            GMGInfo<Prec_> & _info;

        public:
            unsigned long current_level;
            GMGState(GMGInfo<Prec_> & info) :
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
                coarse_smoother_functor();
            }

    };

    template<typename Tag_, typename Mode_>
    struct GMG
    {
        public:
            template<typename Prec_>
            static void value(GMGInfo<Prec_> & info)
            {
                GMGState<Prec_> state(info);
                std::list<unsigned long>::iterator cycle_iterator(info.cycle.begin());
                std::list<std::tr1::function<void ()> >::iterator smoother_functors_iterator(info.smoother_functors.begin());
                std::list<std::tr1::function<void ()> >::iterator coarse_solver_functors_iterator(info.coarse_solver_functors.begin());
                std::list<std::tr1::function<void ()> >::iterator transfer_functors_iterator(info.transfer_functors.begin());
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
