/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/product.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_product(DenseVector<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
        {
            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }
            if (a.rows() != result.size())
            {
                throw VectorSizeDoesNotMatch(a.rows(), result.size());
            }

            void * b_cl(b.lock(lm_read_only, Tag_::memory_value));
            void * result_cl(result.lock(lm_write_only, Tag_::memory_value));
            void * Aj_cl(a.Aj().lock(lm_read_only, Tag_::memory_value));
            void * Ax_cl(a.Ax().lock(lm_read_only, Tag_::memory_value));
            void * Arl_cl(a.Arl().lock(lm_read_only, Tag_::memory_value));
            std::string opname("product_smell_dv_");
            opname += typeid(DT_).name();
            opencl::product_smell_dv<DT_>(b_cl, result_cl, Aj_cl, Ax_cl, Arl_cl,
                    a.rows(), a.columns(), a.num_cols_per_row(), a.stride(), a.threads(), tag_to_device<Tag_>(), opname);
            result.unlock(lm_write_only);
            b.unlock(lm_read_only);
            a.Aj().unlock(lm_read_only);
            a.Ax().unlock(lm_read_only);
            a.Arl().unlock(lm_read_only);
        }
    }
}

using namespace honei;

template <typename DT_>
DenseVector<DT_> & Product<tags::OpenCL::CPU>::value(DenseVector<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<DT_> with DenseVectorContinuousBase<DT_> (OpenCL CPU):");
    PROFILER_START("Product tags::OpenCL::CPU");
    opencl::common_product<tags::OpenCL::CPU>(result, a, b);
    PROFILER_STOP("Product tags::OpenCL::CPU");
    return result;
}
template DenseVector<float> & Product<tags::OpenCL::CPU>::value(DenseVector<float> &, const SparseMatrixELL<float> &, const DenseVector<float> &);
template DenseVector<double> & Product<tags::OpenCL::CPU>::value(DenseVector<double> &, const SparseMatrixELL<double> &, const DenseVector<double> &);

template <typename DT_>
DenseVector<DT_> & Product<tags::OpenCL::GPU>::value(DenseVector<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
{
    CONTEXT("When multiplying SparseMatrixELL<DT_> with DenseVectorContinuousBase<DT_> (OpenCL GPU):");
    PROFILER_START("Product tags::OpenCL::GPU");
    opencl::common_product<tags::OpenCL::GPU>(result, a, b);
    PROFILER_STOP("Product tags::OpenCL::GPU");
    return result;
}
template DenseVector<float> & Product<tags::OpenCL::GPU>::value(DenseVector<float> &, const SparseMatrixELL<float> &, const DenseVector<float> &);
template DenseVector<double> & Product<tags::OpenCL::GPU>::value(DenseVector<double> &, const SparseMatrixELL<double> &, const DenseVector<double> &);
