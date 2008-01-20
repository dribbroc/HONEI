/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/assertion.hh>
#include <honei/libutil/memory_backend_cell.hh>
#include <honei/libutil/spe_manager.hh>
#include <honei/libutil/stringify.hh>
#include <unittest/unittest.hh>

#include <cmath>
#include <map>

#include <honei/libspe2.h>

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




            CellBackend::instance()->download(id_a, spe.id(), c, 32 * sizeof(float));
            CellBackend::instance()->download(id_b, spe.id(), c, 32 * sizeof(float));

            CellBackend::instance()->swap(id_a, id_b);
            CellBackend::instance()->download(id_a, spe.id(), c, 32 * sizeof(float));
            CellBackend::instance()->download(id_b, spe.id(), c, 32 * sizeof(float));

            CellBackend::instance()->free(id_a, spe.id());
            CellBackend::instance()->free(id_b, spe.id());
        }
} cell_backend_function_test;


