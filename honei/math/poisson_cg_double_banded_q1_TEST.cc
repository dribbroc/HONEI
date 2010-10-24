/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/math/jacobi.hh>
#include <honei/math/conjugate_gradients.hh>
#include <honei/math/iterative_refinement.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <honei/swe/post_processing.hh>
//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class PoissonTestCGBandedFloatQ1:
    public BaseTest
{
    public:
        PoissonTestCGBandedFloatQ1(const std::string & tag) :
            BaseTest("Poisson test for itrerative LES solvers , CG (Q1 system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long n(25);
            DenseVector<double> dd_v(n, double(0));
            DenseVector<double> ll_v(n, double(0));
            DenseVector<double> ld_v(n, double(0));
            DenseVector<double> lu_v(n, double(0));
            DenseVector<double> dl_v(n, double(0));
            DenseVector<double> du_v(n, double(0));
            DenseVector<double> ul_v(n, double(0));
            DenseVector<double> ud_v(n, double(0));
            DenseVector<double> uu_v(n, double(0));
            DenseVector<double> b_v(n, double(0));
            DenseVector<double> ana_sol_v(n, double(0));
            DenseVector<double> ref_sol_v(n, double(0));
            long root_n = (long)sqrt(n);
            BandedMatrixQ1<double> A(n,ll_v, ld_v , lu_v, dl_v, dd_v, du_v, ul_v, ud_v, uu_v);

            FillMatrix<Tag_, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(A);
            FillVector<Tag_, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(b_v);
            DenseVector<double> result(n, double(0));
            ConjugateGradients<Tag_, NONE>::value(A, b_v, result, 20ul);
            result.lock(lm_read_only);
            DenseMatrix<double> view(root_n, root_n);
            unsigned long r(0), c(0);
            for(unsigned long index(0) ; index < n ; ++ index)
            {
                view[r][c] = result[index];
                if((index + 1) % root_n == 0)
                {
                    ++r;
                    c=0;
                }
                else
                    ++c;
            }
            //PostProcessing<output_types::GNUPLOT>::value(view, 1, root_n, root_n, 0);
            result.unlock(lm_read_only);
            TEST_CHECK(true);
        }
};
PoissonTestCGBandedFloatQ1<tags::CPU, double> poisson_test_cg_banded_double("double");
/*#ifdef HONEI_SSE
PoissonTestCGBandedFloatQ1<tags::CPU::SSE, double> poisson_test_cg_banded_double_sse("SSE double");
#endif
#ifdef HONEI_CUDA
PoissonTestCGBandedFloatQ1<tags::GPU::CUDA, double> poisson_test_cg_banded_double_cuda("double");
#endif
#ifdef HONEI_CELL
PoissonTestCGBandedFloatQ1<tags::Cell, double> poisson_test_cg_banded_double_cell("Cell double");
#endif*/
