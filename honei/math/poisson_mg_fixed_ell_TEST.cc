/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
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

#define SOLVER_VERBOSE_L2
#define SOLVER_VERBOSE

#include <honei/math/multigrid.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class PoissonTestMGSparseELL:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f;

    public:
        PoissonTestMGSparseELL(const std::string & tag,
                unsigned long size,
                std::string file) :
            BaseTest("Poisson test for itrerative LES solvers , MG (ELLPACK system)<" + tag + ">")
    {
        register_tag(Tag_::name);
        _size = size;
        _res_f = file;
    }
        virtual void run() const
        {
            unsigned long _root_n(_size);
            unsigned long n(_root_n * _root_n);
            MGInfo<DT1_, SparseMatrixELL<DT1_> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;
            DenseVector<unsigned long> mask(8);

            info.macro_border_mask = new DenseVector<unsigned long>(8);
            for(int i(0); i < 8; ++i)
            {
                (*info.macro_border_mask)[i] = 2;
            }
            //set Neumann boundaries:
            (*info.macro_border_mask)[5] =1;

            info.min_level = 1;
            switch(n)
            {
                case 1050625:
                    {
                        info.max_level = 10;
                    }
                    break;
                case 263169:
                    {
                        info.max_level = 9;
                    }
                    break;
                case 66049:
                    {
                        info.max_level = 8;
                    }
                    break;
                case 16641:
                    {
                        info.max_level = 7;
                    }
                    break;
                case 4225:
                    {
                        info.max_level = 6;
                    }
                    break;
                case 1089:
                    {
                        info.max_level = 5;
                    }
                    break;
                case 289:
                    {
                        info.max_level = 4;
                    }
                    break;
                case 81:
                    {
                        info.max_level = 3;
                    }
                    break;
                case 25:
                    {
                        info.max_level = 2;
                    }
                    break;
                case 9:
                    {
                        info.max_level = 1;
                    }
                    break;
                default:
                    throw InternalError("Uknown size!");
                    break;
            }

            info.n_max_iter = 16;
            info.initial_zero = false;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 2;
            info.n_post_smooth = 2;
            info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.tolerance_coarse = 1e-2;
            info.adapt_correction_factor = 1.;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                BandedMatrixQ1<DT1_> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<DT1_> sm(ac_a);
                SparseMatrixELL<DT1_> ac_s(sm);
                info.a.push_back(ac_s);
                // iteration vectors
                DenseVector<DT1_> ac_c(size, DT1_(0));
                info.c.push_back(ac_c);
                DenseVector<DT1_> ac_d(size, DT1_(0));
                info.d.push_back(ac_d);
                DenseVector<DT1_> ac_rhs(size, DT1_(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<DT1_> ac_x(size, DT1_(0));
                info.x.push_back(ac_x);

                info.diags_inverted.push_back(dummy_band.copy());
            }


            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                std::cout << size << std::endl;
                // iteration vectors
                DenseVector<DT1_> ac_c(size, DT1_(0));
                info.c.push_back(ac_c);
                DenseVector<DT1_> ac_d(size, DT1_(0));
                info.d.push_back(ac_d);
                DenseVector<DT1_> ac_x(size, DT1_(0));
                info.x.push_back(ac_x);

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                //info.diags_inverted.push_back(dummy_band.copy());
            }

            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                DenseVector<DT1_> band(N);
                BandedMatrixQ1<DT1_> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                DenseVector<DT1_> current_rhs(N);

                FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_matrix);

                FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_rhs);

                info.rhs.push_back(current_rhs);

                SparseMatrix<DT1_> sm(current_matrix);
                SparseMatrixELL<DT1_> smell(sm);
                info.a.push_back(smell);

                DenseVector<DT1_> scaled_diag_inverted(current_matrix.band(DD).copy());
                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(size==0)
                    size = 9;

                DenseVector<DT1_> null(size , DT1_(0));
                info.x[i] = null.copy();
            }

            DenseVector<DT1_> result(n, DT1_(0));
            DenseVector<DT1_> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<DT1_> system(info.a[info.max_level]);
            Multigrid<Tag_, Tag_, JAC, CYCLE::V, FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
        }
};
PoissonTestMGSparseELL<tags::CPU, float> poisson_test_mg_banded_float("float", 33ul, "1089.bin");
PoissonTestMGSparseELL<tags::CPU, double> poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#ifdef HONEI_SSE
PoissonTestMGSparseELL<tags::CPU::SSE, float> sse_poisson_test_mg_banded_float("float", 33ul, "1089.bin");
PoissonTestMGSparseELL<tags::CPU::MultiCore::SSE, float> mc_sse_poisson_test_mg_banded_float("float", 33ul, "1089.bin");
PoissonTestMGSparseELL<tags::CPU::SSE, double> sse_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
PoissonTestMGSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#endif
#ifdef HONEI_CUDA
PoissonTestMGSparseELL<tags::GPU::CUDA, float> cuda_poisson_test_mg_banded_float("float", 33ul, "1089.bin");
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGSparseELL<tags::GPU::CUDA, double> cuda_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#endif
#endif
