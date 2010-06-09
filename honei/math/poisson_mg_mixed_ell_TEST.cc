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
#include <honei/math/vector_io.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/transposition.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename ITag_, typename OTag_, typename DT1_>
class PoissonTestMGSparseELLMixed:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f;

    public:
        PoissonTestMGSparseELLMixed(const std::string & tag,
                unsigned long size,
                std::string file) :
            BaseTest("Poisson test for itrerative LES solvers , MG (ELLPACK system)<" + tag + ">")
    {
        register_tag(OTag_::name);
        _size = size;
        _res_f = file;
    }
        virtual void run() const
        {
            unsigned long _root_n(_size);
            unsigned long n(_root_n * _root_n);
            MGInfo<float, SparseMatrixELL<float> > info;
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

            info.n_max_iter = 4;
            info.initial_zero = false;
            info.tolerance = 1e-8;
            info.convergence_check = false;

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

                DenseVector<float> dummy_band(size, float(0));
                BandedMatrixQ1<float> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<float> sm(ac_a);
                SparseMatrixELL<float> ac_s(sm);
                info.a.push_back(ac_s);
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_rhs(size, float(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                info.diags_inverted.push_back(dummy_band.copy());
            }


            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                std::cout << size << std::endl;
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                DenseVector<float> dummy_band(size, float(0));
                //info.diags_inverted.push_back(dummy_band.copy());
            }

            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                DenseVector<float> band(N);
                BandedMatrixQ1<float> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                DenseVector<float> current_rhs(N);

                FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_matrix);

                FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_rhs);

                info.rhs.push_back(current_rhs);

                SparseMatrix<float> sm(current_matrix);
                SparseMatrixELL<float> smell(sm);
                info.a.push_back(smell);

                DenseVector<float> scaled_diag_inverted(current_matrix.band(DD).copy());
                ElementInverse<OTag_>::value(scaled_diag_inverted);
                Scale<OTag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
                info.x[i] = null.copy();
            }

            DenseVector<DT1_> null(info.rhs[info.max_level].size() , DT1_(0));
            BandedMatrixQ1<DT1_> A(info.rhs[info.max_level].size() , null.copy(), null.copy() , null.copy(), null.copy() , null.copy(), null.copy(), null.copy(), null.copy(), null.copy());
            FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(A);
            SparseMatrix<DT1_> sm(A);
            SparseMatrixELL<DT1_> system(sm);
            DenseVector<DT1_> RHS( info.rhs[info.max_level].size(), DT1_(0.));
            FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(RHS);
            DenseVector<DT1_> result(n, DT1_(0));
            Multigrid<ITag_, OTag_, NONE, JAC, CYCLE::V, MIXED >::value(system, RHS, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
        }
};
/*PoissonTestMGSparseELLMixed<tags::CPU, tags::CPU, double> poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#ifdef HONEI_SSE
PoissonTestMGSparseELLMixed<tags::CPU::SSE, tags::CPU::SSE, double> sse_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
PoissonTestMGSparseELLMixed<tags::CPU::MultiCore::SSE, tags::CPU::MultiCore::SSE, double> mcsse_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_SSE
PoissonTestMGSparseELLMixed<tags::GPU::CUDA, tags::CPU::SSE, double> cuda_poisson_test_mg_banded_double_2("double", 33ul, "1089.bin");
PoissonTestMGSparseELLMixed<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> mcsse_cuda_poisson_test_mg_banded_double_2("double", 33ul, "1089.bin");
#endif
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGSparseELLMixed<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_test_mg_banded_double("double", 33ul, "1089.bin");
#endif
#endif*/

template <typename ITag_, typename OTag_, typename DT1_>
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
        register_tag(OTag_::name);
        _size = size;
    }
        virtual void run() const
        {
            unsigned long _root_n(_size);
            unsigned long n(_root_n * _root_n);
            MGInfo<float, SparseMatrixELL<float> > info;
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
                    throw InternalError("Unknown size!");
                    break;
            }

            info.min_level = 2;
            info.n_max_iter = 1000;
            info.initial_zero = true;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            info.n_max_iter_coarse = ((unsigned long)sqrt((float)(pow((float)2 , (float)info.max_level) + 1)*(pow((float)2 , (float)info.max_level) + 1)));
            info.tolerance_coarse = std::numeric_limits<double>::epsilon();
            info.adapt_correction_factor = 1.;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1)));
                if(i == 0)
                    size = 9;

                DenseVector<float> dummy_band(size, float(0));
                BandedMatrixQ1<float> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<float> sm(ac_a);
                SparseMatrixELL<float> ac_s(sm);
                info.a.push_back(ac_s);
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_rhs(size, float(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                info.diags_inverted.push_back(dummy_band.copy());
            }


            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size = (unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1));
                std::cout << size << std::endl;
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                DenseVector<float> dummy_band(size, float(0));
                //info.diags_inverted.push_back(dummy_band.copy());
            }

            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N = (unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1));
                DenseVector<float> current_rhs(N);
                std::string A_file(HONEI_SOURCEDIR);
                A_file += "/honei/math/testdata/poisson/";
                A_file += "poisson_A_";
                A_file += stringify(i);
                A_file += ".ell";
                SparseMatrixELL<float> smell(MatrixIO<io_formats::ELL, SparseMatrixELL<double> >::read_matrix(A_file, float(0)));

                std::string rhs_file(HONEI_SOURCEDIR);
                rhs_file += "/honei/math/testdata/poisson/";
                rhs_file += "poisson_rhs";
                if(i == info.max_level)
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, float(0));

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<float> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<ITag_>::value(scaled_diag_inverted);
                Scale<ITag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<float> prol(1, 1);
                        SparseMatrixELL<float> prolmat(prol);
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
                        SparseMatrixELL<float> prolmat(MatrixIO<io_formats::ELL, SparseMatrixELL<double> >::read_matrix(prol_file, float(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<float> prol(prolmat);
                        SparseMatrix<float> res(prol.columns(), prol.rows());
                        Transposition<ITag_>::value(prol, res);
                        SparseMatrixELL<float> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((float)2, (float)i) + 1) * ((unsigned long)pow((float)2, (float)i) + 1)));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
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
            DenseVector<float> result(n, float(0));
            DenseVector<float> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<float> system(info.a[info.max_level]);
            Multigrid<ITag_, OTag_, methods::PROLMAT, JAC, CYCLE::V, MIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<float>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
            std::string sol_file(HONEI_SOURCEDIR);
            sol_file += "/honei/math/testdata/poisson/";
            sol_file += "poisson_sol";
            DenseVector<float> ref_result(VectorIO<io_formats::EXP>::read_vector(sol_file, float(0)));

            for(unsigned long i(0) ; i < ref_result.size() ; ++i)
            {
                std::cout << result[i] << " " << ref_result[i] << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], std::numeric_limits<double>::epsilon()*1e14);
            }
        }
};
/*PoissonTestMGSparseELLProlMat<tags::CPU, tags::CPU, double> poisson_test_mg_sparse_prolmat_double_1("double", 65ul);
#ifdef HONEI_SSE
PoissonTestMGSparseELLProlMat<tags::CPU::SSE,tags::CPU::SSE, double> sse_poisson_test_mg_sparse_prolmat_double_4("double", 65ul);
PoissonTestMGSparseELLProlMat<tags::CPU::MultiCore::SSE, tags::CPU::MultiCore::SSE, double> mc_sse_poisson_test_mg_sparse_prolmat_double_5("double", 65ul);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_test_mg_sparse_prolmat_double_8("double", 65ul);
#endif
#endif
#if defined HONEI_CUDA && defined HONEI_SSE
#if defined HONEI_CUDA_DOUBLE
PoissonTestMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::SSE, double> cuda_poisson_test_mg_sparse_prolmat_double_11("double", 65ul);
PoissonTestMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cuda_poisson_test_mg_sparse_prolmat_double_12("double", 65ul);
#endif
#endif
*/

template <typename ITag_, typename OTag_, typename DT1_>
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
            BaseTest("Poisson advanced mixedprec test, ProlMat, MG (ELLPACK system)<" + tag + ">" + " Level= " + stringify(level) + " Sorting= " + stringify(sorting))
    {
        register_tag(OTag_::name);
        _size = level;
        _sorting = sorting;
        _file_base = file_base;
        _nc = nc;
    }
        virtual void run() const
        {
            unsigned long n(_level_to_size(_size, _nc));
            MGInfo<DT1_, SparseMatrixELL<DT1_> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;
            DenseVector<unsigned long> mask(8);

            /*info.macro_border_mask = new DenseVector<unsigned long>(8);
            for(int i(0); i < 8; ++i)
            {
                (*info.macro_border_mask)[i] = 2;
            }*/
            //set Neumann boundaries:
            //(*info.macro_border_mask)[5] =1;


            info.min_level = 2;
            info.max_level = _size;
            info.n_max_iter = 10000;
            info.initial_zero = false;
            info.tolerance = 1e-2;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            //info.n_max_iter_coarse = ((unsigned long)sqrt((DT1_)(pow((DT1_)2 , (DT1_)info.max_level) + 1)*(pow((DT1_)2 , (DT1_)info.max_level) + 1)));
            info.n_max_iter_coarse = 10000;
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
            std::cout << "File:" << file_base << std::endl;
            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N(_level_to_size(i, _nc));
                DenseVector<DT1_> current_rhs(N);
                std::string A_file(file_base);
                A_file += "A_";
                A_file += stringify(i);
                A_file += ".ell";
                SparseMatrixELL<DT1_> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

                std::string rhs_file(file_base);
                rhs_file += "rhs_" + stringify(_size);
                if(i == info.max_level)
                    VectorIO<io_formats::EXP>::read_vector(rhs_file, current_rhs);

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<DT1_> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<OTag_>::value(scaled_diag_inverted);
                Scale<OTag_>::value(scaled_diag_inverted, 0.7);

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
                        SparseMatrixELL<DT1_> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, DT1_(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<DT1_> prol(prolmat);
                        SparseMatrix<DT1_> res(prol.columns(), prol.rows());
                        Transposition<OTag_>::value(prol, res);
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
            DenseVector<DT1_> result(n, DT1_(0));
            VectorIO<io_formats::EXP>::read_vector(init_file, result);

            DenseVector<DT1_> rhs(info.rhs[info.max_level]);
            SparseMatrixELL<DT1_> system(info.a[info.max_level]);
            Multigrid<ITag_, OTag_, methods::PROLMAT, JAC, CYCLE::V, MIXED>::value(system, rhs, result, (unsigned long)11, std::numeric_limits<DT1_>::epsilon(), info);
            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            //std::cout<< result <<endl;
            DenseVector<DT1_> ref_result(rhs.size());
            std::string sol_file(file_base);
            sol_file += "sol_" + stringify(_size);
            VectorIO<io_formats::EXP>::read_vector(sol_file, ref_result);

            for(unsigned long i(0) ; i < ref_result.size() ; ++i)
            {
                std::cout << result[i] << " " << ref_result[i] << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], std::numeric_limits<DT1_>::epsilon()*1e12);
            }
        }
};
#ifdef HONEI_SSE
//  PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::SSE, tags::CPU::SSE, double> sse_poisson_advanced_test_mg_sparse_prolmat_double("double", 4ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
//  PoissonAdvancedTestMGSparseELLProlMat<tags::CPU::MultiCore::SSE, tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_test_mg_sparse_prolmat_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
//  PoissonAdvancedTestMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_test_mg_sparse_prolmat_double("double", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif
