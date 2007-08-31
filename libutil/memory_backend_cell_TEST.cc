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
#include <map>

#include <libspe2.h>

using namespace honei;
using namespace tests;

class CellBackendFunctionTest :
    public BaseTest
{
    public:
        CellBackendFunctionTest() :
            BaseTest("cell_backend_function_test")
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

            CellBackend::instance()->free(id_a, spe.id());
            CellBackend::instance()->free(id_b, spe.id());
        }
} cell_backend_function_test;

class CellBackendAlignmentTest :
    public QuickTest
{
    public:
        CellBackendAlignmentTest() :
            QuickTest("cell_backend_alignment_test")
        {
        }

        virtual void run() const
        {
            char aligned_data1[32] __attribute__ ((aligned (16)));
            char * aligned_data2(new char[32]);
            SPE spe(*SPEManager::instance()->begin());

            MemoryId id1(10), id2(11);

            bool has_thrown(false);
            try
            {
                CellBackend::instance()->upload(id1, spe.id(), aligned_data1, 32);
                CellBackend::instance()->upload(id2, spe.id(), aligned_data2, 32);
            }
            catch (Exception & e)
            {
                has_thrown = true;
                throw;
            }
            TEST_CHECK(! has_thrown);

            char * misaligned_data1(aligned_data1 + 3);
            char * misaligned_data2(aligned_data2 + 7);

            MemoryId id3(12), id4(13);
            TEST_CHECK_THROWS(CellBackend::instance()->upload(id3, spe.id(), misaligned_data1, 16),
                    MisalignmentError);
            TEST_CHECK_THROWS(CellBackend::instance()->upload(id4, spe.id(), misaligned_data2, 16),
                    MisalignmentError);

        }
} cell_backend_alignment_test;
