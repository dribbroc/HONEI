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

//#define SOLVER_VERBOSE_L2
//#define SOLVER_VERBOSE

#include <honei/math/multigrid.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <honei/math/transposition.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/math/endian_swap.hh>
#include <stdio.h>
#include <stdlib.h>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/spai.hh>

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
                unsigned long size) :
            BaseTest("Poisson test for itrerative LES solvers , MG (ELLPACK system)<" + tag + ">")
    {
        register_tag(Tag_::name);
        _size = size;
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
            Multigrid<Tag_, Tag_, methods::NONE, methods::JAC, methods::CYCLE::V, methods::FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
        }
};
/*PoissonTestMGSparseELL<tags::CPU, double> poisson_test_mg_sparse_double("double", 33ul);
#ifdef HONEI_SSE
PoissonTestMGSparseELL<tags::CPU::SSE, double> sse_poisson_test_mg_sparse_double("double", 33ul);
PoissonTestMGSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_poisson_test_mg_sparse_double("double", 33ul);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGSparseELL<tags::GPU::CUDA, double> cuda_poisson_test_mg_sparse_double("double", 33ul);
#endif
#endif*/

template <typename Tag_, typename DT1_>
class PoissonTestMGSparseELLProlMat:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f;

    public:
        PoissonTestMGSparseELLProlMat(const std::string & tag,
                unsigned long size) :
            BaseTest("Poisson test for itrerative LES solvers, ProlMat, MG (ELLPACK system)<" + tag + ">")
    {
        register_tag(Tag_::name);
        _size = size;
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
            //(*info.macro_border_mask)[5] =1;

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

            info.min_level = 1;
            info.n_max_iter = 1000;
            info.initial_zero = true;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.tolerance_coarse = std::numeric_limits<double>::epsilon();
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
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
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
                DenseVector<DT1_> current_rhs(N);
                std::string A_file(HONEI_SOURCEDIR);
                A_file += "/honei/math/testdata/poisson/";
                A_file += "poisson_A_";
                A_file += stringify(i);
                A_file += ".ell";
                SparseMatrixELL<DT1_> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

                std::string rhs_file(HONEI_SOURCEDIR);
                rhs_file += "/honei/math/testdata/poisson/";
                rhs_file += "poisson_rhs";
                if(i == info.max_level)
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0));

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<DT1_> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<DT1_> prol(1, 1);
                        SparseMatrixELL<DT1_> prolmat(prol);
                        info.prolmats.push_back(prolmat);
                        info.resmats.push_back(prolmat);
                    }
                    else
                    {
                        std::string prol_file(HONEI_SOURCEDIR);
                        prol_file += "/honei/math/testdata/poisson/";
                        prol_file += "poisson_prol_";
                        prol_file += stringify(i);
                        prol_file += ".ell";
                        SparseMatrixELL<DT1_> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, DT1_(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<DT1_> prol(prolmat);
                        SparseMatrix<DT1_> res(prol.columns(), prol.rows());
                        Transposition<Tag_>::value(prol, res);
                        SparseMatrixELL<DT1_> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
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

            /*for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                std::cout << "VECSIZE d " << info.d.at(i).size() << std::endl;
                std::cout << "VECSIZE rhs " << info.rhs.at(i).size() << std::endl;
                std::cout <<"SYSTEM"<<std::endl;
                std::cout <<"ROW: " << info.a.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.a.at(i).columns() << std::endl;
                std::cout <<"RESTRICTION"<<std::endl;
                std::cout <<"ROW: " << info.resmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.resmats.at(i).columns() << std::endl;
                std::cout <<"PROLONGATION"<<std::endl;
                std::cout <<"ROW: " << info.prolmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.prolmats.at(i).columns() << std::endl;

            }*/
            DenseVector<DT1_> result(n, DT1_(0));
            DenseVector<DT1_> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<DT1_> system(info.a[info.max_level]);
            Multigrid<Tag_, Tag_, methods::PROLMAT, methods::JAC, methods::CYCLE::V, methods::FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
            std::string sol_file(HONEI_SOURCEDIR);
            sol_file += "/honei/math/testdata/poisson/";
            sol_file += "poisson_sol";
            DenseVector<DT1_> ref_result(VectorIO<io_formats::EXP>::read_vector(sol_file, DT1_(0)));

            for(unsigned long i(0) ; i < ref_result.size() ; ++i)
            {
                std::cout << result[i] << " " << ref_result[i] << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], std::numeric_limits<DT1_>::epsilon()*1e11);
            }
        }
};
/*PoissonTestMGSparseELLProlMat<tags::CPU, double> poisson_test_mg_sparse_prolmat_double("double", 65ul);
#ifdef HONEI_SSE
  PoissonTestMGSparseELLProlMat<tags::CPU::SSE, double> sse_poisson_test_mg_sparse_prolmat_double("double", 65ul);
  PoissonTestMGSparseELLProlMat<tags::CPU::MultiCore::SSE, double> mc_sse_poisson_test_mg_sparse_prolmat_double("double", 65ul);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGSparseELLProlMat<tags::GPU::CUDA, double> cuda_poisson_test_mg_sparse_prolmat_double("double", 65ul);
#endif
#endif*/

template <typename Tag_, typename DT1_>
class PoissonAdvancedTestMGSparseELLProlMat:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;
        unsigned _nc;

        static unsigned long _level_to_size(unsigned long level, unsigned elem_type) // 0 = q1 , 1 = q1t , 2 = q2
        {
            switch(level)
            {
                case 10:
                    {
                        if(elem_type == 0)
                            return 2101248;
                    }
                case 9:
                    {
                        if(elem_type == 0)
                            return 526336;
                        else if(elem_type == 1)
                            return 1050624;
                        else
                            return 2101248;
                    }
                case 8:
                    {
                        if(elem_type == 0)
                            return 132096;
                        else if(elem_type == 1)
                            return 263168;
                        else
                            return 526336;
                    }
                case 7:
                    {
                        if(elem_type == 0)
                            return 33280;
                        else if(elem_type == 1)
                            return 66048;
                        else
                            return 132096;
                    }
                case 6:
                    {
                        if(elem_type == 0)
                            return 8448;
                        else if(elem_type == 1)
                            return 16640;
                        else
                            return 33280;
                    }
                case 5:
                    {
                        if(elem_type == 0)
                            return 2176;
                        else if(elem_type == 1)
                            return 4224;
                        else
                            return 8448;
                    }
                case 4:
                    {
                        if(elem_type == 0)
                            return 576;
                        else if(elem_type == 1)
                            return 1088;
                        else
                            return 2176;
                    }
                case 3:
                    {
                        if(elem_type == 0)
                            return 160;
                        else if(elem_type == 1)
                            return 288;
                        else
                            return 576;
                    }
                case 2:
                    {
                        if(elem_type == 0)
                            return 48;
                        else if(elem_type == 1)
                            return 80;
                        else
                            return 160;
                    }
                case 1:
                    {
                        if(elem_type == 0)
                            return 16;
                        else if(elem_type == 1)
                            return 24;
                        else
                            return 48;
                    }
                default:
                    return 1;
            }
        }

    public:
        PoissonAdvancedTestMGSparseELLProlMat(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            BaseTest("Poisson advanced test for itrerative LES solvers, ProlMat, MG (ELLPACK system)<" + tag + ">" + " Level= " + stringify(level) + " Sorting= " + stringify(sorting))
    {
        register_tag(Tag_::name);
        _size = level;
        _sorting = sorting;
        _file_base = file_base;
        _nc = nc;
    }
        virtual void run() const
        {
            //unsigned long n(_level_to_size(_size, _nc));
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
            //(*info.macro_border_mask)[5] =1;


            info.min_level = 1;
            info.max_level = _size;
            info.n_max_iter = 30;
            info.initial_zero = false;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            //info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.n_max_iter_coarse = 1000;
            info.tolerance_coarse = 1e-8;
            info.adapt_correction_factor = 1.;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(i == 0)
                    size = 9;

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                BandedMatrixQ1<DT1_> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<DT1_> sm(ac_a);
                SparseMatrixELL<DT1_> ac_s(sm);
                info.a.push_back(ac_s);
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
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
                unsigned long size(_level_to_size(i, _nc));
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

            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";
            std::cout << "Filebase:" << file_base << std::endl;
            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N(_level_to_size(i, _nc));
                DenseVector<DT1_> current_rhs(N);
                std::string A_file(file_base);
                A_file += "A_";
                A_file += stringify(i);
                A_file += ".ell";
                std::cout << "File:" << A_file << std::endl;
                SparseMatrixELL<DT1_> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

                std::string rhs_file(file_base);
                rhs_file += "rhs_" + stringify(_size);
                std::cout << "File:" << rhs_file << std::endl;
                if(i == info.max_level)
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0));

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<DT1_> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<DT1_> prol(1, 1);
                        SparseMatrixELL<DT1_> prolmat(prol);
                        info.prolmats.push_back(prolmat);
                        info.resmats.push_back(prolmat);
                    }
                    else
                    {
                        std::string prol_file(file_base);
                        prol_file += "prol_";
                        prol_file += stringify(i);
                        prol_file += ".ell";
                        std::cout << "File:" << prol_file << std::endl;
                        SparseMatrixELL<DT1_> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, DT1_(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<DT1_> prol(prolmat);
                        SparseMatrix<DT1_> res(prol.columns(), prol.rows());
                        Transposition<Tag_>::value(prol, res);
                        SparseMatrixELL<DT1_> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(size==0)
                    size = 9;

                DenseVector<DT1_> null(size , DT1_(0));
                info.x[i] = null.copy();
            }

            /*for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                std::cout << "VECSIZE d " << info.d.at(i).size() << std::endl;
                std::cout << "VECSIZE rhs " << info.rhs.at(i).size() << std::endl;
                std::cout <<"SYSTEM"<<std::endl;
                std::cout <<"ROW: " << info.a.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.a.at(i).columns() << std::endl;
                std::cout <<"RESTRICTION"<<std::endl;
                std::cout <<"ROW: " << info.resmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.resmats.at(i).columns() << std::endl;
                std::cout <<"PROLONGATION"<<std::endl;
                std::cout <<"ROW: " << info.prolmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.prolmats.at(i).columns() << std::endl;

            }*/
            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            std::cout << "File:" << init_file << std::endl;
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            DenseVector<DT1_> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<DT1_> system(info.a[info.max_level]);
            Multigrid<Tag_, Tag_, methods::PROLMAT, methods::JAC, methods::CYCLE::V, methods::FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
            std::string sol_file(file_base);
            sol_file += "sol_" + stringify(_size);
            std::cout << "File:" << sol_file << std::endl;
            DenseVector<DT1_> ref_result(VectorIO<io_formats::EXP>::read_vector(sol_file, DT1_(0)));

            for(unsigned long i(0) ; i < ref_result.size() ; ++i)
            {
                //std::cout << result[i] << " " << ref_result[i] << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], std::numeric_limits<DT1_>::epsilon()*1e12);
            }
        }
};
#ifdef HONEI_SSE
PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::SSE, double> sse_poisson_advanced_test_mg_sparse_prolmat_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);

PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_test_mg_sparse_prolmat_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);

PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::SSE, double> sse_poisson_advanced_test_mg_sparse_prolmat_double_q2("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_test_mg_sparse_prolmat_double_q2("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE

PoissonAdvancedTestMGSparseELLProlMat<tags::GPU::CUDA, double> cuda_poisson_advanced_test_mg_sparse_prolmat_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);

PoissonAdvancedTestMGSparseELLProlMat<tags::GPU::CUDA, double> cuda_poisson_advanced_test_mg_sparse_prolmat_double_q2("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

#endif
#endif

template <typename Tag_, typename DT1_>
class PoissonAdvancedTestMGSparseELLProlMatSpai:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;
        unsigned _nc;

        static unsigned long _level_to_size(unsigned long level, unsigned elem_type) // 0 = q1 , 1 = q1t , 2 = q2
        {
            switch(level)
            {
                case 10:
                    {
                        if(elem_type == 0)
                            return 2101248;
                    }
                case 9:
                    {
                        if(elem_type == 0)
                            return 526336;
                        else if(elem_type == 1)
                            return 1050624;
                        else
                            return 2101248;
                    }
                case 8:
                    {
                        if(elem_type == 0)
                            return 132096;
                        else if(elem_type == 1)
                            return 263168;
                        else
                            return 526336;
                    }
                case 7:
                    {
                        if(elem_type == 0)
                            return 33280;
                        else if(elem_type == 1)
                            return 66048;
                        else
                            return 132096;
                    }
                case 6:
                    {
                        if(elem_type == 0)
                            return 8448;
                        else if(elem_type == 1)
                            return 16640;
                        else
                            return 33280;
                    }
                case 5:
                    {
                        if(elem_type == 0)
                            return 2176;
                        else if(elem_type == 1)
                            return 4224;
                        else
                            return 8448;
                    }
                case 4:
                    {
                        if(elem_type == 0)
                            return 576;
                        else if(elem_type == 1)
                            return 1088;
                        else
                            return 2176;
                    }
                case 3:
                    {
                        if(elem_type == 0)
                            return 160;
                        else if(elem_type == 1)
                            return 288;
                        else
                            return 576;
                    }
                case 2:
                    {
                        if(elem_type == 0)
                            return 48;
                        else if(elem_type == 1)
                            return 80;
                        else
                            return 160;
                    }
                case 1:
                    {
                        if(elem_type == 0)
                            return 16;
                        else if(elem_type == 1)
                            return 24;
                        else
                            return 48;
                    }
                default:
                    return 1;
            }
        }

    public:
        PoissonAdvancedTestMGSparseELLProlMatSpai(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            BaseTest("Poisson advanced test for itrerative LES solvers, ProlMat, MG, Spai (ELLPACK system)<" + tag + ">" + " Level= " + stringify(level) + " Sorting= " + stringify(sorting))
    {
        register_tag(Tag_::name);
        _size = level;
        _sorting = sorting;
        _file_base = file_base;
        _nc = nc;
    }
        virtual void run() const
        {
            //unsigned long n(_level_to_size(_size, _nc));
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
            //(*info.macro_border_mask)[5] =1;


            info.min_level = 1;
            info.max_level = _size;
            info.n_max_iter = 30;
            info.initial_zero = false;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            //info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.n_max_iter_coarse = 1000;
            info.tolerance_coarse = 1e-8;
            info.adapt_correction_factor = 1.;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(i == 0)
                    size = 9;

                DenseVector<DT1_> dummy_band(size, DT1_(0));
                BandedMatrixQ1<DT1_> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<DT1_> sm(ac_a);
                SparseMatrixELL<DT1_> ac_s(sm);
                info.a.push_back(ac_s);
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
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
                info.spais.push_back(ac_s.copy());
            }

            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
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

            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";
            std::cout << "Filebase:" << file_base << std::endl;
            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N(_level_to_size(i, _nc));
                DenseVector<DT1_> current_rhs(N);
                std::string A_file(file_base);
                A_file += "A_";
                A_file += stringify(i);
                A_file += ".ell";
                std::cout << "File:" << A_file << std::endl;
                SparseMatrixELL<DT1_> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

                if(i == info.max_level)
                {
                    std::string rhs_file(file_base);
                    rhs_file += "rhs_" + stringify(_size);
                    std::cout << "File:" << rhs_file << std::endl;
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0));
                }

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<DT1_> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<Tag_>::value(scaled_diag_inverted);
                Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                SparseMatrix<DT1_> spai(smell);
                double stellen(log10(smell.used_elements()));
                double spai_eps(stellen / 10.);
                std::cout<<"Calc spai for lvl "<<i<<" with eps="<<spai_eps<<std::endl;
                SparseMatrix<DT1_> sm_m(SPAI::value(spai, spai_eps, 35, 35));
                SparseMatrixELL<DT1_> spai_m(sm_m);
                info.spais.push_back(spai_m.copy());

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<DT1_> prol(1, 1);
                        SparseMatrixELL<DT1_> prolmat(prol);
                        info.prolmats.push_back(prolmat);
                        info.resmats.push_back(prolmat);
                    }
                    else
                    {
                        std::string prol_file(file_base);
                        prol_file += "prol_";
                        prol_file += stringify(i);
                        prol_file += ".ell";
                        std::cout << "File:" << prol_file << std::endl;
                        SparseMatrixELL<DT1_> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, DT1_(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<DT1_> prol(prolmat);
                        SparseMatrix<DT1_> res(prol.columns(), prol.rows());
                        Transposition<Tag_>::value(prol, res);
                        SparseMatrixELL<DT1_> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(size==0)
                    size = 9;

                DenseVector<DT1_> null(size , DT1_(0));
                info.x[i] = null.copy();
            }

            /*for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                std::cout << "VECSIZE d " << info.d.at(i).size() << std::endl;
                std::cout << "VECSIZE rhs " << info.rhs.at(i).size() << std::endl;
                std::cout <<"SYSTEM"<<std::endl;
                std::cout <<"ROW: " << info.a.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.a.at(i).columns() << std::endl;
                std::cout <<"RESTRICTION"<<std::endl;
                std::cout <<"ROW: " << info.resmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.resmats.at(i).columns() << std::endl;
                std::cout <<"PROLONGATION"<<std::endl;
                std::cout <<"ROW: " << info.prolmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.prolmats.at(i).columns() << std::endl;

            }*/
            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            std::cout << "File:" << init_file << std::endl;
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            DenseVector<DT1_> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<DT1_> system(info.a[info.max_level]);
            Multigrid<Tag_, Tag_, methods::PROLMAT, methods::SPAI, methods::CYCLE::V, methods::FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
            std::string sol_file(file_base);
            sol_file += "sol_" + stringify(_size);
            std::cout << "File:" << sol_file << std::endl;
            DenseVector<DT1_> ref_result(VectorIO<io_formats::EXP>::read_vector(sol_file, DT1_(0)));

            for(unsigned long i(0) ; i < ref_result.size() ; ++i)
            {
                //std::cout << result[i] << " " << ref_result[i] << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], std::numeric_limits<DT1_>::epsilon()*1e12);
            }
        }
};
#ifdef HONEI_SSE
PoissonAdvancedTestMGSparseELLProlMatSpai<tags::CPU::SSE, double> sse_poisson_advanced_test_mg_sparse_prolmat_spai_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedTestMGSparseELLProlMatSpai<tags::CPU::SSE, double> sse_poisson_advanced_test_mg_sparse_prolmat_spai_double_2("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
