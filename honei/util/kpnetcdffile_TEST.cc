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

#include <honei/util/netcdf_datatypes.hh>
#include <honei/util/kpnetcdffile.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace honei;
using namespace tests;

class KPNetCDFFileTest :
    public QuickTest
{
    public:
        KPNetCDFFileTest() :
            QuickTest("kpnetcdffile_test")
        {
        }

        virtual void run() const
        {
            std::string path(HONEI_SOURCEDIR);
            path += "/honei/lbm/sample.nc";
            KPInitialConditions cond;
            shared_ptr<Field> b;
            shared_ptr<Field> u1;
            shared_ptr<Field> u2;
            shared_ptr<Field> u3;
            KPNetCDFFile file(path, cond, b, u1, u2, u3);
            unsigned long nx, ny;
            TEST_CHECK(cond.isValid());
            nx = cond.getNx();
            ny = cond.getNy();
            std::cout<<"Gridsize: "<<nx<<" x "<<ny<<std::endl;
            std::cout<<u1->data[5]<<std::endl;
            shared_ptr<TimeStep> ts(new TimeStep(nx, ny));
            file.readTimeStep(ts, 0);
            std::cout<<ts->U[0]->data[5]<<std::endl;
            file.readTimeStep(ts, 30);
            std::cout<<ts->U[0]->data[5]<<std::endl;
        }
} kpnetcdffile_test;
