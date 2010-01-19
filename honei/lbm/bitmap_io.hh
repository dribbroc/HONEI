/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LBM_GUARD_BITMAP_IO_HH
#define LBM_GUARD_BITMAP_IO_HH 1

#include <honei/util/pnm_image.hh>
#include <honei/la/dense_matrix.hh>

using namespace honei;
using namespace util;

namespace io_formats
{
    class BMP;
    class PGM;
    class PPM;
}

template<typename FileType_>
class BitmapIO
{
};

template<>
class BitmapIO<io_formats::PGM>
{
    public:
        template<typename DT_>
        static DenseMatrix<DT_> read_scalar_field(std::string source, DT_ scale, int channel=1)
        {
            std::cout << "Reading PGM image from file " << source << ":" << std::endl;
            std::tr1::shared_ptr<PGMImage> image(PGMImage::read(source.c_str()));

            unsigned long height((unsigned long)image->getHeight());
            unsigned long width((unsigned long)image->getWidth());
            std::cout << "Lines         : " << height << std::endl;
            std::cout << "Pixel per line: " << width << std::endl;

            std::cout << "Reading data..." << std::endl;
            float * data(image->getGrayData());
            DenseMatrix<DT_> result(height, width);
            for(unsigned long i(0) ; i < height ; ++i)
                for(unsigned long j(0) ; j < width ; ++j)
                    result[i][j] = (DT_)data[i * width + j];

            std::cout << "...finished." << std::endl;

            return result;
        }
};

template<>
class BitmapIO<io_formats::PPM>
{
    public:
        template<typename DT_>
        static DenseMatrix<DT_> read_scalar_field(std::string source, DT_ scale, int channel=1)
        {
            std::cout << "Reading PGM image from file " << source << ":" << std::endl;
            std::tr1::shared_ptr<PPMImage> image(PPMImage::read(source.c_str()));

            unsigned long height((unsigned long)image->getHeight());
            unsigned long width((unsigned long)image->getWidth());
            std::cout << "Lines         : " << height << std::endl;
            std::cout << "Pixel per line: " << width << std::endl;

            std::cout << "Reading data..." << std::endl;
            float * data;
            switch(channel)
            {
                case 1:
                    {
                        data = image->getRedData();
                    }
                    break;

                case 2:
                    {
                        data = image->getGreenData();
                    }

                case 3:
                    {
                        data = image->getBlueData();
                    }
            }
            DenseMatrix<DT_> result(height, width);
            for(unsigned long i(0) ; i < height ; ++i)
                for(unsigned long j(0) ; j < width ; ++j)
                    result[i][j] = (DT_)data[i * width + j];

            std::cout << "...finished." << std::endl;

            return result;
        }
};
#endif
