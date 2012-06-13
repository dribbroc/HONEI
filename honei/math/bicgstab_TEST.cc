/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. LibMath is free software;
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

#include <honei/math/bicgstab.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <iomanip>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class BiCGStabSolverTestSparseELLPrecon:
    public BaseTest
{
    public:
        BiCGStabSolverTestSparseELLPrecon(const std::string & tag):
            BaseTest("BiCGStabSolver solver test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {

            std::string base("/home/user/dribbroc/sandkasten/jonas/");
            std::string filename(base);
            filename += "A.txt";
            SparseMatrix<DT1_> ssmatrix2(MatrixIO<io_formats::MTX>::read_matrix(filename, DT1_(0)));
            SparseMatrixELL<DT1_> smatrix2(ssmatrix2);

            std::string filename_2(base);
            filename_2 += "b.txt";
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));

            DenseVector<DT1_> diag_inverted(smatrix2.rows(), DT1_(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = DT1_(1);//DT1_(0.7)/smatrix2(i, i);
            }

            std::string filename_4(base);
            filename_4 += "x.txt";
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(filename_4, DT1_(0)));
            unsigned long used_iters;
            BiCGStabSolver<Tag_, methods::VAR>::value(smatrix2, diag_inverted, rhs, result, 700ul, used_iters, DT1_(1e-12));

            result.lock(lm_read_only);
            //std::cout << result << std::endl;
            result.unlock(lm_read_only);

            std::cout << "Used iters: " << used_iters << std::endl;
        }
};
BiCGStabSolverTestSparseELLPrecon<tags::CPU, double> cg_precon_test_double_sparse_ell("double JAC");
