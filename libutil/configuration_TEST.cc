/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <libutil/configuration.hh>
#include <unittest/unittest.hh>

#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace honei;
using namespace tests;

class ConfigurationTest :
    public QuickTest
{
    public:
        ConfigurationTest() :
            QuickTest("configuration_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::string filename("/tmp/.honeirc");

                std::remove(filename.c_str());

                {
                    CONTEXT("When creating configuration file:");
                    std::fstream file(filename.c_str(), std::ios_base::out);

                    file << "number-of-cores=1" << std::endl;
                    file << "key-with-trailing-whitespaces \t     =2" << std::endl;
                    file << "value-with-leading-whitespaces= \t    3" << std::endl;
                    file << "key-and-value-with-whitespaces \t=\t  4" << std::endl;

                    file.close();
                }

                ::setenv("HONEI_CONFIG", filename.c_str(), 1);

                {
                    CONTEXT("When reading configuration data:");
                    // This may read the configuration twice, but it is needed.
                    Configuration::instance()->reread();

                    TEST_CHECK_EQUAL(Configuration::instance()->get_value("number-of-cores", 0), 1);
                    TEST_CHECK_EQUAL(Configuration::instance()->get_value("key-with-trailing-whitespaces", 1), 2);
                    TEST_CHECK_EQUAL(Configuration::instance()->get_value("value-with-leading-whitespaces", 2), 3);
                    TEST_CHECK_EQUAL(Configuration::instance()->get_value("key-and-value-with-whitespaces", 3), 4);
                }

                std::remove(filename.c_str());
                TEST_CHECK(true);
            }
            catch (Exception & e)
            {
                std::cout << e.backtrace("\n ...") << std::endl;
                TEST_CHECK(false);
            }
            catch (std::exception & e)
            {
                std::cout << e.what() << std::endl;
                TEST_CHECK(false);
            }
        }
} configuration_test;
