/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <benchmark/benchmark.hh>
#include <honei/math/matrix_io.hh>

#include <cstdlib>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class ELLMatrixIOBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
        std::string _file;

    public:
        ELLMatrixIOBench(const std::string & id, unsigned long size, int count, std::string file) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
            _file = file;
        }


        virtual void run()
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _file;
            BENCHMARK(
                    SparseMatrixELL<DataType_> smatrix = MatrixIO<io_formats::ELL>::read_matrix(filename, DataType_(1));
                    );
            evaluate();
        }
};
ELLMatrixIOBench<tags::CPU, double>  ELLMIOBench0 ("ELL: 5pt_10x10.ell",  1, 1, "5pt_10x10.ell");
//ELLMatrixIOBench<tags::CPU, double>  ELLMIOBench1 ("ELL: rrze1.ell",  1, 1, "rrze1.ell");

template <typename Tag_, typename DataType_>
class MTXMatrixIOBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
        std::string _file;

    public:
        MTXMatrixIOBench(const std::string & id, unsigned long size, int count, std::string file) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
            _file = file;
        }



        virtual void run()
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _file;
            {
                BENCHMARK(
                        for (unsigned long i(0) ; i < 1 ; ++i)
                        {
                        unsigned long non_zeros(MatrixIO<io_formats::MTX>::get_non_zeros(filename));
                        unsigned long rows;
                        unsigned long columns;
                        unsigned long ax;
                        unsigned long bx;
                        DenseVector<unsigned long> r(non_zeros);
                        DenseVector<unsigned long> c(non_zeros);
                        DenseVector<double> data(non_zeros);

                        MatrixIO<io_formats::MTX>::read_matrix(filename, r, c, data);
                        MatrixIO<io_formats::MTX>::get_sizes(filename, rows, columns, ax, bx);
                        SparseMatrix<double> tsmatrix(rows, columns, r, c, data);
                        SparseMatrixELL<double> smatrix(tsmatrix);
                        }
                        );
            }
            evaluate();
        }
};
MTXMatrixIOBench<tags::CPU, double>  MTXMIOBench0 ("MTX: 5pt_10x10.mtx",  1, 1, "5pt_10x10.mtx");
//MTXMatrixIOBench<tags::CPU, double>  MTXMIOBench1 ("MTX: rrze1.mtx",  1, 1, "rrze1.mtx");
