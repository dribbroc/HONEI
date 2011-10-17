/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
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


#include <honei/math/fft.hh>
#include <honei/math/ri.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <math.h>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT_>
class SmootherTestSparse:
    public BaseTest
{
    public:
        SmootherTestSparse(const std::string & tag) :
            BaseTest("Smoother test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            //Read in data:
            std::string dir(HONEI_SOURCEDIR);
            std::string file (dir + "/honei/math/testdata/poisson_advanced2/q1t_sort_3/");
            file += "A_6";
            file += ".ell";
            SparseMatrixELL<DT_> A(MatrixIO<io_formats::ELL>::read_matrix(file, DT_(0)));

            std::string file2 (dir + "/honei/math/testdata/poisson_advanced2/q1t_sort_3/");
            file2 += "rhs_6";
            DenseVector<DT_> rhs(VectorIO<io_formats::EXP>::read_vector(file2, DT_(0)));

            /*DenseVector<DT_> diag_inverted(A.rows(), DT_(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = DT_(0.7)/A(i, i);
            }*/
            std::string file5 (dir + "/honei/math/testdata/poisson_advanced2/q1t_sort_3/");
            file5 += "A_6_spai_grote";
            file5 += ".ell";
            SparseMatrixELL<DT_> diag_inverted(MatrixIO<io_formats::ELL>::read_matrix(file5, DT_(0)));

            std::string file3 (dir + "/honei/math/testdata/poisson_advanced2/q1t_sort_3/");
            file3 += "init_6";
            DenseVector<DT_> result(VectorIO<io_formats::EXP>::read_vector(file2, DT_(0)));

            DenseVector<DT_> x(rhs.size());
            Defect<Tag_>::value(x, rhs, A, result);
            DenseVector<DT_> pre_y = FFT::forward(x);

            DenseVector<DT_> temp_0(rhs.size());
            DenseVector<DT_> temp_1(rhs.size());
            RISmoother<Tag_>::value(A, diag_inverted, rhs, result, temp_0, temp_1, 100ul);

            Defect<Tag_>::value(x, rhs, A, result);
            DenseVector<DT_> post_y = FFT::forward(x);

            DenseVector<DT_> quotient(post_y.size());
            for (unsigned long i(0) ; i < post_y.size() ; ++i)
            {
                quotient[i] = (pre_y[i] == DT_(0)) || fabs(post_y[i] / pre_y[i]) > 1 ? DT_(0) : fabs(post_y[i] / pre_y[i]);
            }


            /*std::string out_filename (dir + "/honei/math/testdata/");
            out_filename += "fft.dat";
            std::ofstream hfile;
            hfile.open(out_filename.c_str(), fstream::trunc | fstream::out);

            for (unsigned long i(0) ; i < quotient.size() ; ++i)
            {
                std::string record = stringify(i) + " " + stringify(quotient[i]) + "\n";
                hfile << record;
            }
            hfile.close();*/

            // http://de.wikipedia.org/wiki/Kleinste-Quadrate-Sch%C3%A4tzer#Spezialfall_einer_einfachen_linearen_Ausgleichsgeraden
            DT_ mean_x(0);
            DT_ mean_y(0);
            for (unsigned long i(0) ; i < quotient.size() ; ++i)
            {
                mean_x += i;
                mean_y += quotient[i];
            }
            mean_x /= quotient.size();
            mean_y /= quotient.size();

            DenseVector<DT_> dx(quotient.size());
            DenseVector<DT_> dy(quotient.size());
            for (unsigned long i(0) ; i < quotient.size() ; ++i)
            {
                dx[i] = i - mean_x;
                dy[i] = quotient[i] - mean_y;
            }

            DT_ sum_above(0);
            DT_ sum_beneath(0);
            for (unsigned long i(0) ; i < quotient.size() ; ++i)
            {
                sum_above += dx[i] * dy[i];
                sum_beneath += dx[i] * dx[i];
            }
            DT_ alpha(sum_above / sum_beneath);
            DT_ alpha_0(mean_y - (alpha * mean_x));

            std::cout<<"alpha: "<<alpha<<" alpha_0: "<<alpha_0<<std::endl;

            std::string out_filename (dir + "/honei/math/testdata/");
            out_filename += "fft.dat";
            std::ofstream hfile;
            hfile.open(out_filename.c_str(), fstream::trunc | fstream::out);

            for (unsigned long i(0) ; i < quotient.size() ; ++i)
            {
                std::string record = stringify(i) + " " + stringify(alpha_0 + (DT_(i) * alpha)) + "\n";
                hfile << record;
            }
            hfile.close();
        }
};
SmootherTestSparse<tags::CPU, double> smoother_test_double("double");
