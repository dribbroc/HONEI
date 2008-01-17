/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 , 2008 Markus Geveler <apryde@gmx.de>
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

#include <libmath/jacobi_kernel.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;


template <typename Tag_, typename DT1_>
class JacobiKernelTestBanded:
    public BaseTest
{
    public:
        JacobiKernelTestBanded(const std::string & tag) :
            BaseTest("Jacobi (SINGLE STEPPING) solver test (Banded system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> diagd(4, DT1_(1));
            DenseVector<DT1_> bp1(4, DT1_(1));
            DenseVector<DT1_> bp2(4, DT1_(1));
            DenseVector<DT1_> bm1(4, DT1_(1));
            DenseVector<DT1_> bm2(4, DT1_(1));
            DenseVector<DT1_> bp3(4, DT1_(0));
            DenseVector<DT1_> bm3(4, DT1_(0));

            diagd[0] = DT1_(7);
            diagd[1] = DT1_(6);
            diagd[2] = DT1_(5);
            diagd[3] = DT1_(1);

            bp1[0] = DT1_(-2);
            bp1[1] = DT1_(2);
            bp1[2] = DT1_(0);
            bp1[3] = DT1_(0);

            bp2[0] = DT1_(0);
            bp2[1] = DT1_(0);
            bp2[2] = DT1_(0);
            bp2[3] = DT1_(0);

            bm1[0] = DT1_(0);
            bm1[1] = DT1_(-2);
            bm1[2] = DT1_(2);
            bm1[3] = DT1_(0);

            bm2[0] = DT1_(0);
            bm2[1] = DT1_(0);
            bm2[2] = DT1_(0);
            bm2[3] = DT1_(0);

            BandedMatrix<DT1_> A(4);
            A.insert_band(0, diagd);
            A.insert_band(1, bp1);
            A.insert_band(2, bp2);
            A.insert_band(-1, bm1);
            A.insert_band(-2, bm2);
            A.insert_band(3, bp3);
            A.insert_band(-3, bm3);

            DenseVector<DT1_> b(4, DT1_(1));
            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);
            b[3] = DT1_(1);

            std::cout<<"A:" << A <<endl;
            std::cout<<"b:" << b <<endl;

            DenseVector<DT1_> x_analytical(4, DT1_(1));
            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            x_analytical[3] = DT1_(1.);

            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            DenseVector<DT1_> x(b.size(), DT1_(0));
            DenseVector<DT1_> x_last(x.copy());

            DT1_ norm_x_last = DT1_(0);
            DT1_ norm_x = DT1_(1);
            DenseVector<DT1_> diag(b.size(), DT1_(0));

            DenseVector<DT1_> diag_inverted(b.size(), DT1_(0));

            BandedMatrix<DT1_> difference(A.copy());
            ///Create Diagonal, invert, compute difference on the fly.
            for(unsigned long i =0; i < diag.size(); ++i)
            {
                diag[i] = A.band(0)[i];
                if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                {
                    diag_inverted[i] = DT1_(1) / diag[i];
                }
                else
                {
                    diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                }
            }

            DenseVector<DT1_> zeros(b.size(), DT1_(0));
            difference.insert_band(0, zeros);
            //Scale<tags::CPU>::value(difference, DT1_(-1));

            DT1_ konv_rad = std::numeric_limits<DT1_>::epsilon();
            while(fabs(norm_x - norm_x_last) > konv_rad)
            {

                JacobiKernel<Tag_>::value(b, x, diag_inverted, difference);
                norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                x_last = x.copy();
            }
            cout << "x: " << x << endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, norm_x , 10e-06);

        }
};

template <typename Tag_, typename DT1_>
class JacobiKernelTestBandedSSE:
    public BaseTest
{
    public:
        JacobiKernelTestBandedSSE(const std::string & tag) :
            BaseTest("Jacobi (SINGLE STEPPING) solver test (Banded system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> diagd(4, DT1_(1));
            DenseVector<DT1_> bp1(4, DT1_(1));
            DenseVector<DT1_> bp2(4, DT1_(1));
            DenseVector<DT1_> bm1(4, DT1_(1));
            DenseVector<DT1_> bm2(4, DT1_(1));
            DenseVector<DT1_> bp3(4, DT1_(0));
            DenseVector<DT1_> bm3(4, DT1_(0));

            diagd[0] = DT1_(7);
            diagd[1] = DT1_(6);
            diagd[2] = DT1_(5);
            diagd[3] = DT1_(1);

            bp1[0] = DT1_(-2);
            bp1[1] = DT1_(2);
            bp1[2] = DT1_(0);
            bp1[3] = DT1_(0);

            bp2[0] = DT1_(0);
            bp2[1] = DT1_(0);
            bp2[2] = DT1_(0);
            bp2[3] = DT1_(0);

            bm1[0] = DT1_(0);
            bm1[1] = DT1_(-2);
            bm1[2] = DT1_(2);
            bm1[3] = DT1_(0);

            bm2[0] = DT1_(0);
            bm2[1] = DT1_(0);
            bm2[2] = DT1_(0);
            bm2[3] = DT1_(0);

            BandedMatrix<DT1_> A(4);
            A.insert_band(0, diagd);
            A.insert_band(1, bp1);
            A.insert_band(2, bp2);
            A.insert_band(-1, bm1);
            A.insert_band(-2, bm2);
            A.insert_band(3, bp3);
            A.insert_band(-3, bm3);

            DenseVector<DT1_> b(4, DT1_(1));
            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);
            b[3] = DT1_(1);

            std::cout<<"A:" << A <<endl;
            std::cout<<"b:" << b <<endl;

            DenseVector<DT1_> x_analytical(4, DT1_(1));
            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            x_analytical[3] = DT1_(1.);

            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            DenseVector<DT1_> x(b.size(), DT1_(0));
            DenseVector<DT1_> x_last(x.copy());

            DT1_ norm_x_last = DT1_(0);
            DT1_ norm_x = DT1_(1);
            DenseVector<DT1_> diag(b.size(), DT1_(0));

            DenseVector<DT1_> diag_inverted(b.size(), DT1_(0));

            BandedMatrix<DT1_> difference(A.copy());
            ///Create Diagonal, invert, compute difference on the fly.
            for(unsigned long i =0; i < diag.size(); ++i)
            {
                diag[i] = A.band(0)[i];
                if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                {
                    diag_inverted[i] = DT1_(1) / diag[i];
                }
                else
                {
                    diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                }
            }

            DenseVector<DT1_> zeros(b.size(), DT1_(0));
            difference.insert_band(0, zeros);
            Scale<tags::CPU>::value(difference, DT1_(-1));

            DT1_ konv_rad = std::numeric_limits<DT1_>::epsilon();
            while(fabs(norm_x - norm_x_last) > konv_rad)
            {

                JacobiKernel<Tag_>::value(b, x, diag_inverted, difference);
                norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                x_last = x.copy();
            }
            cout << "x: " << x << endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, norm_x , 10e-06);

        }
};

JacobiKernelTestBanded<tags::CPU, float> jacobi_kernel_test_float_banded("float");
JacobiKernelTestBanded<tags::CPU, double> jacobi_kernel_test_double_banded("double");
#ifdef HONEI_SSE
JacobiKernelTestBandedSSE<tags::CPU::SSE, float> jacobi_kernel_test_float_banded_sse("SSE float");
JacobiKernelTestBandedSSE<tags::CPU::SSE, double> jacobi_kernel_test_double_banded_sse("SSE double");

#endif

