/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBMATH_GUARD_QUADRATURE_HH
#define LIBMATH_GUARD_QUADRATURE_HH 1

#include <libutil/tags.hh>
#include <functional>

namespace honei
{

    namespace tags
    {
        struct Trapezoid;
    }

    template <typename Tag_, typename F_> class Quadrature;

    /**
     * Create and return suitable Quadrature functor for a given function
     * and stepsize.
     *
     * \param function The function that shall be integrated.
     * \param step_size The step size that will be used for integration.
     */
    template <typename Tag_, typename F_> inline Quadrature<Tag_, F_> make_quadrature(const F_ & function,
            const typename F_::argument_type & step_size)
    {
        return Quadrature<Tag_, F_>(function, step_size);
    }

    /**
     * Quadrature<tags::Trapezoid> uses the trapezoid method to integrate a
     * unary function in a given range with a given step size.
     */
    template <typename F_> class Quadrature<tags::Trapezoid, F_> :
        public std::binary_function<typename F_::argument_type,
               typename F_::argument_type,
               typename F_::result_type>
    {
        private:
            /// Our function that will be integrated.
            const F_ _function;

            /// Our step size.
            const typename F_::argument_type _step_size;

            /// Constructor.
            Quadrature(const F_ & function, const typename F_::argument_type 
                    & step_size) :
                _function(function),
                _step_size(step_size)
            {
            }

        public:
            friend Quadrature<tags::Trapezoid, F_> make_quadrature<tags::Trapezoid, F_>(const F_ &,
                    const typename F_::argument_type &);

            /**
             * Integrate a previously specified function in the range [left, right]
             * using the trapezoid method.
             *
             * \param left Left border of the quadrature interval.
             * \param right Right border of the quadrature interval.
             */
            typename F_::result_type operator() (const typename F_::argument_type & left,
                    const typename F_::argument_type & right)
            {
                typename F_::result_type result(0);

                for (typename F_::argument_type a(left), b(left + _step_size) ;
                        a < right ; a = b, b += _step_size)
                {
                    result += _function(a) + _function(b);
                }

                result *= _step_size;
                result /= 2;

                return result;
            }
    };
}

#endif
