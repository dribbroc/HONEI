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

#include <honei/util/bitmap.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>

using namespace honei;
class BitmapIO
{
    public:
        template<typename DT_>
        static DenseMatrix<DT_> read_scalar_field(std::string source)
        {
            CBitmap bitmap(source.c_str());
            unsigned long width(bitmap.GetWidth());
            unsigned long height(bitmap.GetHeight());
            unsigned long bits(bitmap.GetBitCount());

            std::cout << "Loading data from Bitmap..." << std::endl;
            std::cout << "Width    = " << width << std::endl;
            std::cout << "Height   = " << height << std::endl;
            std::cout << "Bitdepth = " << bits << std::endl;

            DenseMatrix<DT_> result(height, width);
            void* data_v(bitmap.GetBits());
            RGBA* data = (RGBA*)data_v;


            for(unsigned long i(0) ; i < width * height ; ++i)
            {
                std::cout << "R: " << (float)data[i].Red << std::endl;
                std::cout << "G: " << (float)data[i].Green << std::endl;
                std::cout << "B: " << (float)data[i].Blue << std::endl;
                std::cout << "A: " << (float)data[i].Alpha << std::endl;
            }
            return result;
        }
};

#endif
