/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <libutil/assertion.hh>
#include <libutil/memory_backend_cell.hh>
#include <libutil/spe_manager.hh>
#include <libutil/stringify.hh>
#include <unittest/unittest.hh>

#include <cmath>
#include <iostream>
#include <map>

#include <libspe2.h>

using namespace pg512;
using namespace tests;

class CellBackendTest :
    public BaseTest
{
    public:
        CellBackendTest() :
            BaseTest("cell_backend_test")
        {
        }

        virtual void run() const
        {
            float * a(new float[32]),
                  * b(new float[32]),
                  * c(new float[32]);
            SPE spe(*SPEManager::instance()->begin());

            std::fill_n(a, 32, float(1.11));
            std::generate_n(b, 32, ::rand);

            MemoryId id_a(10), id_b(11);
            CellBackend::instance()->upload(id_a, spe.id(), a, 32 * sizeof(float));
            CellBackend::instance()->upload(id_b, spe.id(), b, 32 * sizeof(float));

            const CellBackend::Chunk * chunk_a(CellBackend::instance()->get_chunk(id_a, spe.id())),
                * chunk_b(CellBackend::instance()->get_chunk(id_b, spe.id()));

            TEST_CHECK(32 * sizeof(float) == chunk_a->size);
            TEST_CHECK(32 * sizeof(float) == chunk_b->size);

            char * ls(static_cast<char *>(spe_ls_area_get(spe.context())));
            TEST_CHECK(0 == memcmp(ls + chunk_a->address, a, 32 * sizeof(float)));
            TEST_CHECK(0 == memcmp(ls + chunk_b->address, b, 32 * sizeof(float)));

            CellBackend::instance()->download(id_a, spe.id(), c, 32 * sizeof(float));
            TEST_CHECK(0 == memcmp(c, a, 32 * sizeof(float)));
            CellBackend::instance()->download(id_b, spe.id(), c, 32 * sizeof(float));
            TEST_CHECK(0 == memcmp(c, b, 32 * sizeof(float)));

            CellBackend::instance()->swap(id_a, id_b);
            CellBackend::instance()->download(id_a, spe.id(), c, 32 * sizeof(float));
            TEST_CHECK(0 == memcmp(c, b, 32 * sizeof(float)));
            CellBackend::instance()->download(id_b, spe.id(), c, 32 * sizeof(float));
            TEST_CHECK(0 == memcmp(c, a, 32 * sizeof(float)));
        }
} memory_manager_test;
