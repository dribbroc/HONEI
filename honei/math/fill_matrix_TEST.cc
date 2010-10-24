/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. Math is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * Math is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/math/fill_matrix.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

#include <limits>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class FillMatrixTestDirichlet0:
    public BaseTest
{
    private:
        unsigned long _size;
    public:
        FillMatrixTestDirichlet0(const std::string & tag, unsigned long size) :
            BaseTest("FillMatrixTestDirichlet0 <" + tag + "> with size: " + stringify(size))
        {
            this->_size = size;
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            //load reference data:
            int n;
            FILE* file;

            double* dd;

            double* ll;
            double* ld;
            double* lu;
            double* dl;
            double* du;
            double* ul;
            double* ud;
            double* uu;

            double* b;
            double* ana_sol;
            double* ref_sol;

            std::string file_name(HONEI_SOURCEDIR);
            file_name += "/honei/math/testdata/" + stringify(_size) + ".bin";
            file = fopen(file_name.c_str(), "rb");
            if (1 != (int)fread(&n, sizeof(int), 1, file))
                throw InternalError("IO Error!");
#ifdef HONEI_CELL
            unsigned char b1, b2, b3, b4;
            b1 = n & 255;
            b2 = ( n >> 8 ) & 255;
            b3 = ( n>>16 ) & 255;
            b4 = ( n>>24 ) & 255;
            n = ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
#endif
            dd = new double[n];
            ll = new double[n];
            ld = new double[n];
            lu = new double[n];
            dl = new double[n];
            du = new double[n];
            ul = new double[n];
            ud = new double[n];
            uu = new double[n];
            b = new double[n];
            ana_sol = new double[n];
            ref_sol = new double[n];

            if (n != (int)fread(dd, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ll, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ld, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(lu, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(dl, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(du, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ul, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ud, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(uu, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(b,  sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ana_sol, sizeof(double), n, file))
                throw InternalError("IO Error!");
            if (n != (int)fread(ref_sol, sizeof(double), n, file))
                throw InternalError("IO Error!");
            fclose(file);
#ifdef HONEI_CELL
            for(int i(0); i < n; ++i)
            {
                dd[i] = DoubleSwap(dd[i]);
                ll[i] = DoubleSwap(ll[i]);
                ld[i] = DoubleSwap(ld[i]);
                lu[i] = DoubleSwap(lu[i]);
                dl[i] = DoubleSwap(dl[i]);
                du[i] = DoubleSwap(du[i]);
                ul[i] = DoubleSwap(ul[i]);
                ud[i] = DoubleSwap(ud[i]);
                uu[i] = DoubleSwap(uu[i]);
                b[i] = DoubleSwap(b[i]);
                ana_sol[i] = DoubleSwap(ana_sol[i]);
                ref_sol[i] = DoubleSwap(ref_sol[i]);
            }
#endif
            DenseVector<DT1_> dd_v(n, DT1_(0));
            DenseVector<DT1_> ll_v(n, DT1_(0));
            DenseVector<DT1_> ld_v(n, DT1_(0));
            DenseVector<DT1_> lu_v(n, DT1_(0));
            DenseVector<DT1_> dl_v(n, DT1_(0));
            DenseVector<DT1_> du_v(n, DT1_(0));
            DenseVector<DT1_> ul_v(n, DT1_(0));
            DenseVector<DT1_> ud_v(n, DT1_(0));
            DenseVector<DT1_> uu_v(n, DT1_(0));
            DenseVector<DT1_> b_v(n, DT1_(0));
            DenseVector<DT1_> ana_sol_v(n, DT1_(0));
            DenseVector<DT1_> ref_sol_v(n, DT1_(0));
            for(int i = 0; i < n; ++i)
            {
                dd_v[i] = (DT1_)dd[i];
                ll_v[i] = (DT1_)ll[i];
                ld_v[i] = (DT1_)ld[i];
                lu_v[i] = (DT1_)lu[i];
                dl_v[i] = (DT1_)dl[i];
                du_v[i] = (DT1_)du[i];
                ul_v[i] = (DT1_)ul[i];
                ud_v[i] = (DT1_)ud[i];
                uu_v[i] = (DT1_)uu[i];
                b_v[i] = (DT1_)b[i];
                ana_sol_v[i] = (DT1_)ana_sol[i];
                ref_sol_v[i] = (DT1_)ref_sol[i];
            }

            BandedMatrixQ1<DT1_> A_ref(n ,ll_v, ld_v , lu_v, dl_v, dd_v, du_v, ul_v, ud_v, uu_v);
            DenseVector<DT1_> null(_size, DT1_(0));

            BandedMatrixQ1<DT1_> result(_size, null.copy(), null.copy(),null.copy(),null.copy(),null.copy(),null.copy(),null.copy(),null.copy(),null.copy());
            FillMatrix<Tag_, applications::POISSON, boundary_types::DIRICHLET::DIRICHLET_0>::value(result);
            result.lock(lm_read_only);
            for(unsigned long i(0) ; i < result.band(DD).size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result.band(DD)[i], A_ref.band(DD)[i], std::numeric_limits<DT1_>::epsilon()*12.);
            }
            result.unlock(lm_read_only);
        }
};
FillMatrixTestDirichlet0<tags::CPU, float> fill_matrix_float_2("float", 25ul);
FillMatrixTestDirichlet0<tags::CPU, float> fill_matrix_float_3("float", 81ul);
FillMatrixTestDirichlet0<tags::CPU, float> fill_matrix_float_4("float", 289ul);
FillMatrixTestDirichlet0<tags::CPU, float> fill_matrix_float_5("float", 1089ul);
FillMatrixTestDirichlet0<tags::CPU, float> fill_matrix_float_6("float", 4225ul);
FillMatrixTestDirichlet0<tags::CPU, double> fill_matrix_double_2("double", 25ul);
FillMatrixTestDirichlet0<tags::CPU, double> fill_matrix_double_3("double", 81ul);
FillMatrixTestDirichlet0<tags::CPU, double> fill_matrix_double_4("double", 289ul);
FillMatrixTestDirichlet0<tags::CPU, double> fill_matrix_double_5("double", 1089ul);
FillMatrixTestDirichlet0<tags::CPU, double> fill_matrix_double_6("double", 4225ul);
