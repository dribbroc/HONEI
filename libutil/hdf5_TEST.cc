/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <libutil/hdf5.hh>
#include <unittest/unittest.hh>

#include <cstdio>
#include <iostream>
#include <sys/stat.h>

using namespace honei;
using namespace tests;

class HDF5CreationTest :
    public QuickTest
{
    public:
        HDF5CreationTest() :
            QuickTest("hdf5_file_creation_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::string filename("hdf5_TEST.h5");

                std::remove(filename.c_str());

                {
                    CONTEXT("When creating hdf5 file:");
                    HDF5File file(filename, 0);
                }

                struct stat buffer;
                TEST_CHECK_EQUAL(0, ::stat(filename.c_str(), &buffer));

                std::remove(filename.c_str());
                TEST_CHECK(true);
            }
            catch (Exception & e)
            {
                std::cout << e.backtrace("\n ...") << std::endl;
                TEST_CHECK(false);
            }
            catch (...)
            {
                TEST_CHECK(false);
            }
        }
} hdf5_creation_test;

class HDF5WriteTest :
    public QuickTest
{
    public:
        HDF5WriteTest() :
            QuickTest("hdf5_write_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::string filename("hdf5_TEST.h5");

                std::remove(filename.c_str());

                {
                    CONTEXT("When creating hdf5 file and writing to it:");
                    HDF5File file(filename, 0);

                    HDF5SimpleDataSpace data_space(1);
                    data_space[10];

                    float data[10];
                    std::fill_n(data, 10, 1.0f);
                    HDF5DataSet<float> data_set(file, "elements", data_space);
                    data_set << data;
                }

                struct stat buffer;
                TEST_CHECK_EQUAL(0, ::stat(filename.c_str(), &buffer));

                std::remove(filename.c_str());
                TEST_CHECK(true);
            }
            catch (Exception & e)
            {
                std::cout << e.backtrace("\n ...") << std::endl;
                TEST_CHECK(false);
            }
            catch (...)
            {
                TEST_CHECK(false);
            }
        }
} hdf5_write_test;
