/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#include <honei/math/bi_conjugate_gradients_stabilised.hh>
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
class SequenceLGSTEST:
    public BaseTest
{
    private:
        std::string _base_f;
    public:
        SequenceLGSTEST(const std::string & tag,
                std::string base_file) :
            BaseTest("Sequence of LGS test" + tag + ">")
        {
            register_tag(Tag_::name);
            _base_f = base_file;
        }

        virtual void run() const
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/testdata/";
            filebase += _base_f;
            filebase += "/";

            for (unsigned long ts(0) ; ts < 5 ; ++ts)
            {
                std::string matrix_f(filebase);
                matrix_f += "A_";
                matrix_f += stringify(ts);
                matrix_f += ".ell";
                SparseMatrixELL<DT1_> smatrix2(MatrixIO<io_formats::ELL>::read_matrix(matrix_f, DT1_(0)));

                std::string rhs_f(filebase);
                rhs_f += "rhs_";
                rhs_f += stringify(ts);
                DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(rhs_f, DT1_(0)));

                std::string precon_f(filebase);
                precon_f += "A_";
                precon_f += stringify(ts);
                precon_f += "_spai.ell";
                SparseMatrixELL<DT1_> smatrix3(MatrixIO<io_formats::ELL>::read_matrix(precon_f, DT1_(0)));

                /*SparseMatrix<DT1_> tsmatrix3(smatrix2.rows(), smatrix2.columns());
                  for(unsigned long i(0) ; i < smatrix2.rows() ; ++i)
                  {
                  tsmatrix3(i , i) = 0.7 * 1/smatrix2(i, i);
                  }
                  SparseMatrixELL<DT1_> smatrix3(tsmatrix3);
                  */

                std::string init_f(filebase);
                init_f += "init_";
                init_f += stringify(ts);
                DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_f, DT1_(0)));

                unsigned long used_iters(0);
                PBiCGStab<Tag_, methods::VAR>::value(smatrix2, rhs, result, smatrix3, 1000ul, used_iters, DT1_(1e-8));
                std::cout<<"timestep " << ts << ", iters: " << used_iters << std::endl;
            }
        }
};
SequenceLGSTEST<tags::CPU, double> pbicgstab_test_double_sparse_ell("double", "sequence_lgs");
