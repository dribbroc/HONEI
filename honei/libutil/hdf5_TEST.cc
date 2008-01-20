/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/hdf5.hh>
#include <unittest/unittest.hh>

#include <cstdio>
#include <iostream>
#include <sys/stat.h>

using namespace honei;
using namespace tests;

class HDF5FileCreationTest :
    public QuickTest
{
    public:
        HDF5FileCreationTest() :
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
                    CONTEXT("When creating HDF5 file:");
                    HDF5File::create(filename, 0);
                }

                struct stat buffer;
                TEST_CHECK_EQUAL(0, ::stat(filename.c_str(), &buffer));

                {
                    CONTEXT("When opening previously created HDF5 file:");
                    HDF5File file(filename);
                }

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
} hdf5_file_creation_test;

class HDF5DataSetTest :
    public QuickTest
{
    public:
        HDF5DataSetTest() :
            QuickTest("hdf5_data_set_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::string filename("hdf5_TEST.h5");

                std::remove(filename.c_str());

                float data[10];
                float data_ref[10];
                std::fill_n(data, 10, 0.0f);
                std::generate_n(data_ref, 10, ::rand);

                {
                    CONTEXT("When creating HDF5 file and writing to it:");
                    HDF5File file(HDF5File::create(filename, 0));

                    HDF5SimpleDataSpace data_space(1);
                    data_space[10];

                    HDF5DataSet<float> data_set(HDF5DataSet<float>::create(file, "elements", data_space));
                    data_set << data_ref;
                }

                struct stat buffer;
                TEST_CHECK_EQUAL(0, ::stat(filename.c_str(), &buffer));

                {
                    CONTEXT("When opening HDF5 file and reading from it:");
                    HDF5File file(filename);

                    HDF5DataSet<float> data_set(file, "elements");
                    data_set >> data;
                }

                for (float * i(data), * i_end(data + 10), * j(data) ; i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<float>::epsilon());
                }

                std::remove(filename.c_str());
                TEST_CHECK(true);
            }
            catch (Exception & e)
            {
                std::cout << e.backtrace("\n ...") << std::endl;
                std::cout << e.what() << std::endl;
                TEST_CHECK(false);
            }
            catch (...)
            {
                std::cout << "Unknown exception" << std::endl;
                TEST_CHECK(false);
            }
        }
} hdf5_data_set_test;

class HDF5RootGroupTest :
    public QuickTest
{
    public:
        HDF5RootGroupTest() :
            QuickTest("hdf5_root_group_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::string filename("hdf5_TEST.h5");

                std::remove(filename.c_str());

                {
                    CONTEXT("When creating HDF5 file:");
                    HDF5File file(HDF5File::create(filename, 0));

                    {
                        CONTEXT("When opening the root group");
                        file.group("/");
                    }
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
} hdf5_root_group_test;
