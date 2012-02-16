/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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


#include <honei/math/sainv.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/la/product.hh>
#include <honei/la/norm.hh>
#include <honei/la/difference.hh>
#include <honei/math/matrix_io.hh>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT_>
class SainvTestSparse:
    public BaseTest
{
    public:
        SainvTestSparse(const std::string & tag) :
            BaseTest("Sainv test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/poisson_advanced/sort_0/A_5.ell";
            SparseMatrixELL<DT_> smell = MatrixIO<io_formats::ELL>::read_matrix(filename, DT_(1));
            SparseMatrix<DT_> sm(smell);
            std::cout<<"Non Zero Elements of A: "<<sm.used_elements()<<std::endl;
            SparseMatrix<DT_> m(SAINV<Tag_>::value(sm));
            std::cout<<"Non Zero Elements of M: "<<m.used_elements()<<std::endl;

            SparseMatrix<DT_> temp(sm.rows(), sm.columns());
            SparseMatrix<DT_> ident(sm.rows(), sm.columns(), 1);
            SparseMatrix<DT_> jac(sm.rows(), sm.columns(), 1);
            for (unsigned long i(0) ; i < ident.rows() ; ++i)
            {
                ident(i, i) = 1;
            }
            for (unsigned long i(0) ; i < ident.rows() ; ++i)
            {
                jac(i, i) = DT_(1) / sm(i, i);
            }
            double min = Norm<vnt_l_one, false, tags::CPU>::value(Difference<tags::CPU>::value(temp, ident, Product<tags::CPU>::value(sm, m)));
            double jacnorm = Norm<vnt_l_one, false, tags::CPU>::value(Difference<tags::CPU>::value(temp, ident, Product<tags::CPU>::value(sm, jac)));
            std::cout<<"SAINV Norm: "<<min<<" Jac Norm: "<<jacnorm<<std::endl;
        }
};
SainvTestSparse<tags::CPU, float> sainv_test_sparse_ell_float("float");
SainvTestSparse<tags::CPU, double> sainv_test_sparse_ell_double("double");
#ifdef HONEI_SSE
SainvTestSparse<tags::CPU::SSE, float> sse_sainv_test_sparse_ell_float("float");
SainvTestSparse<tags::CPU::SSE, double> sse_sainv_test_sparse_ell_double("double");
#endif
