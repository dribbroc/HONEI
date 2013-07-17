/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/lbm/tags.hh>
#include <honei/util/unittest.hh>
#include <honei/lbm/extraction_grid.hh>

#include <limits>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class ExtractionGridTest :
    public TaggedTest<Tag_>
{
    public:
        ExtractionGridTest(const std::string & type) :
            TaggedTest<Tag_>("extraction_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            PackedGridData<lbm_lattice_types::D2Q9, DataType_> data;
            PackedGridInfo<lbm_lattice_types::D2Q9> info;

            data.h = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            info.limits = new DenseVector<unsigned long>(2);
            unsigned long begin(0);
            unsigned long end(data.h->size());
            (*info.limits)[0] = begin;
            (*info.limits)[1] = end;
            data.u = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            data.v = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            data.distribution_x = new DenseVector<DataType_>(9ul, DataType_(2.));
            data.distribution_y = new DenseVector<DataType_>(9ul, DataType_(2.));
            data.f_temp_0 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_1 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_2 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_3 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_4 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_5 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_6 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_7 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_temp_8 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_0 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_1 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_2 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_3 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_4 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_5 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_6 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_7 = new DenseVector<DataType_>(1000, DataType_(1.234));
            data.f_8 = new DenseVector<DataType_>(1000, DataType_(1.234));

            ExtractionGrid<Tag_, lbm_modes::WET>::value(info, data, DataType_(1));


            DenseVector<DataType_> calc_h (data.h->copy());
            DenseVector<DataType_> calc_u (data.u->copy());
            DenseVector<DataType_> calc_v (data.v->copy());

            //for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
            for(unsigned long i(0); i < data.h->size(); ++i)
            {
                (*data.h)[i] = 4711;

                //accumulate
                (*data.h)[i] = (*data.f_0)[i] +
                    (*data.f_1)[i] +
                    (*data.f_2)[i] +
                    (*data.f_3)[i] +
                    (*data.f_4)[i] +
                    (*data.f_5)[i] +
                    (*data.f_6)[i] +
                    (*data.f_7)[i] +
                    (*data.f_8)[i];

                (*data.u)[i] = 4711;
                (*data.u)[i] = ((*data.distribution_x)[0] * (*data.f_0)[i] +
                        (*data.distribution_x)[1] * (*data.f_1)[i] +
                        (*data.distribution_x)[2] * (*data.f_2)[i] +
                        (*data.distribution_x)[3] * (*data.f_3)[i] +
                        (*data.distribution_x)[4] * (*data.f_4)[i] +
                        (*data.distribution_x)[5] * (*data.f_5)[i] +
                        (*data.distribution_x)[6] * (*data.f_6)[i] +
                        (*data.distribution_x)[7] * (*data.f_7)[i] +
                        (*data.distribution_x)[8] * (*data.f_8)[i]) / (*data.h)[i];

                (*data.v)[i] = 4711;
                (*data.v)[i] = ((*data.distribution_y)[0] * (*data.f_0)[i] +
                        (*data.distribution_y)[1] * (*data.f_1)[i] +
                        (*data.distribution_y)[2] * (*data.f_2)[i] +
                        (*data.distribution_y)[3] * (*data.f_3)[i] +
                        (*data.distribution_y)[4] * (*data.f_4)[i] +
                        (*data.distribution_y)[5] * (*data.f_5)[i] +
                        (*data.distribution_y)[6] * (*data.f_6)[i] +
                        (*data.distribution_y)[7] * (*data.f_7)[i] +
                        (*data.distribution_y)[8] * (*data.f_8)[i]) / (*data.h)[i];
            }
            for(unsigned long i(0); i < data.h->size(); ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(calc_h[i], (*data.h)[i], std::numeric_limits<DataType_>::epsilon() * 50.);
                TEST_CHECK_EQUAL_WITHIN_EPS(calc_u[i], (*data.u)[i], std::numeric_limits<DataType_>::epsilon() * 50.);
                TEST_CHECK_EQUAL_WITHIN_EPS(calc_v[i], (*data.v)[i], std::numeric_limits<DataType_>::epsilon() * 50.);
            }

            data.destroy();
            info.destroy();

        }
};
ExtractionGridTest<tags::CPU, float> extraction_grid_test_float("float");
ExtractionGridTest<tags::CPU, double> extraction_grid_test_double("double");
ExtractionGridTest<tags::CPU::Generic, float> generic_extraction_grid_test_float("float");
ExtractionGridTest<tags::CPU::Generic, double> generic_extraction_grid_test_double("double");
#ifdef HONEI_SSE
ExtractionGridTest<tags::CPU::SSE, float> sse_extraction_grid_test_float("float");
ExtractionGridTest<tags::CPU::SSE, double> sse_extraction_grid_test_double("double");
#endif
#ifdef HONEI_CUDA
ExtractionGridTest<tags::GPU::CUDA, float> cuda_extraction_grid_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
ExtractionGridTest<tags::GPU::CUDA, double> cuda_extraction_grid_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
ExtractionGridTest<tags::Cell, float> cell_extraction_grid_test_float("float");
#endif
#ifdef HONEI_OPENCL
ExtractionGridTest<tags::OpenCL::CPU, float> opencl_cpu_extraction_grid_test_float("float");
ExtractionGridTest<tags::OpenCL::CPU, double> opencl_cpu_extraction_grid_test_double("double");
ExtractionGridTest<tags::OpenCL::GPU, float> opencl_gpu_extraction_grid_test_float("float");
ExtractionGridTest<tags::OpenCL::GPU, double> opencl_gpu_extraction_grid_test_double("double");
#endif
