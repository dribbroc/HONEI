/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.tu-dortmund.de>
 *
 * This file is part of the HONEI math C++ library. HONEI is free software;
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

#ifndef MATH_GUARD_PROBLEM_FACTORY_HH
#define MATH_GUARD_PROBLEM_FACTORY_HH 1

#include <vector>
#include <string>
#include <cmath>
#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/math/methods.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/transposition.hh>

#include <iostream>

namespace honei
{
    template<typename DT1_, typename MatrixType_, typename ProlmatType_>
    class FileFactory
    {
        private:
            static std::string _filebase;

        public:
            FileFactory(std::string filebase)
            {
                FileFactory::_filebase = (filebase);
            }

            static void factory(std::vector<MatrixType_ > & a,
                                std::vector<ProlmatType_ > & b,
                                std::vector<ProlmatType_ > & c,
                                std::vector<DenseVector<DT1_> > & d,
                                unsigned long min,
                                unsigned long max)
            {
                /*DenseVector<DT1_> blub(10, 5);
                d.push_back(blub);
                std::cout<<_filebase<<std::endl;*/
                std::string A_file_base(FileFactory::_filebase);
                std::string rhs_file_base(FileFactory::_filebase);
                std::string prol_file_base(FileFactory::_filebase);
                std::string matrix_suffix(".ell");

                A_file_base += "A_";
                rhs_file_base += "rhs";
                prol_file_base += "prol_";

                ///Only store matrices from min to max level
                for(unsigned long i(min) ; i <= max ; ++i)
                {
                    std::string A_file(A_file_base + stringify(i) + matrix_suffix);
                    std::string prol_file(A_file_base + stringify(i) + matrix_suffix);

                    MatrixType_ current_A(MatrixIO<io_formats::ELL, MatrixType_ >::read_matrix(A_file));
                    MatrixType_ current_prol(MatrixIO<io_formats::ELL, MatrixType_ >::read_matrix(prol_file));

                    SparseMatrix<DT1_> prol(current_prol);
                    SparseMatrix<DT1_> res(current_prol.columns(), current_prol.rows());
                    Transposition<tags::CPU>::value(prol, res);
                    SparseMatrixELL<DT1_> current_res(res);

                    a.push_back(current_A);
                    b.push_back(current_prol);
                    c.push_back(current_res);
                }

                DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(rhs_file_base, DT1_(1)));
                d.push_back(rhs);
            }
    };
}

#endif
