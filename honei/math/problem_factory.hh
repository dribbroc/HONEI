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
#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/math/methods.hh>

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
                std::string system_file_base(FileFactory::_filebase);
                std::string rhs_file_base(FileFactory::_filebase);
                std::string prol_file_base(FileFactory::_filebase);

                system_file_base += "_system_";
                rhs_file_base += "_rhs_";
                prol_file_base += "_prol_";

                for(unsigned long i(min) ; i <= max ; ++i)
                {
                    unsigned long n = (unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1));
                    MatrixType_ level_system(MatrixIO<io_formats::ELL, SparseMatrixELL<double> >::read_matrix(A_file, float(0)));
                }
            }
    };
}

#endif
