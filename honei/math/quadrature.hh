/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/util/tags.hh>
#include <functional>
#include <honei/math/interpolation.hh>
#include <honei/libla/dense_matrix.hh>

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

    /**
     * Returns the Volume under a given scalarfield by applying Gaussian Quadrature in 2D.
     *
     */
    template<typename Tag_, typename T_>
    class GaussianQuadrature2D;

    template<typename Tag_>
    class GaussianQuadrature2D<Tag_, tags::Trapezoid>
    {
        private:
            template< typename DT_>
            static DT_ three_point(DenseMatrix<DT_> & s_f, unsigned long i, unsigned long j, unsigned long k, unsigned long l, DT_ delta_x, DT_ delta_y, DT_ right_border, DT_ left_border)
            {
                DT_ h = (right_border - left_border) / s_f.rows();
                DT_ h_1 = (1./2. - 1./(2 *(sqrt(3)))) * h;
                DT_ actual_x = ((k + i) * delta_x) + h_1;
                DT_ actual_y = ((l + j) * delta_y) + h_1;

                DT_ first_point = 5. * Interpolation<Tag_, interpolation_methods::LINEAR>::value(delta_x, delta_y, s_f, actual_x, actual_y);
                actual_x = ((k + i) * delta_x) + h/2.;
                actual_y = ((l + j) * delta_y) + h/2.;
                DT_ second_point = 8. * Interpolation<Tag_, interpolation_methods::LINEAR>::value(delta_x, delta_y, s_f, actual_x, actual_y);

                actual_x = ((k + i) * delta_x) - h_1;
                actual_y = ((l + j) * delta_y) - h_1;
                DT_ third_point = 5. * Interpolation<Tag_, interpolation_methods::LINEAR>::value(delta_x, delta_y, s_f, actual_x, actual_y);

                return first_point + second_point + third_point;
            }

        public:

            /**
             * Returns the Volume under a given scalarfield by applying Gaussian Quadrature in 2D. 
             * \param scalarfield The discrete funktion, that is going to be integrated.
             * \param left_border a in [a,b].
             * \param right_border b in [a,b].
             * \param delta_x The stepsize in x direction.
             * \param delta_y The stepsize in y direction.
             */
            template<typename DT_>
                static inline DT_ value(DenseMatrix<DT_>& scalar_field, DT_ left_border, DT_ right_border, DT_ delta_x, DT_ delta_y)
                {
                    DT_ result(0);
                    for(unsigned long k(0); k < scalar_field.rows(); ++k)
                    {
                        DT_ outermost_sum(0);
                        for(unsigned long l(0); l < scalar_field.rows(); ++l)
                        {
                            DT_ level_sum(0);
                            for(unsigned long i(0); i <= 2; ++i)
                            {
                                DT_ innermost_sum(0);
                                for(unsigned long j(0); j <= 2; ++j)
                                {
                                    innermost_sum += GaussianQuadrature2D<Tag_, tags::Trapezoid>::three_point(scalar_field, i, j, k, l, delta_x, delta_y, right_border, left_border);
                                }
                                innermost_sum *= ((right_border - left_border)/scalar_field.rows())/18;
                                level_sum += innermost_sum;
                            }
                            level_sum *= ((right_border - left_border)/scalar_field.rows())/18;
                            outermost_sum += level_sum;
                        }
                        result += outermost_sum;
                    }
                    return 2 * result;
                }
    };
}

#endif
