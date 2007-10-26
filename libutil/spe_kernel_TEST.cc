/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <cell/cell.hh>
#include <libutil/spe_kernel.hh>
#include <unittest/unittest.hh>

#include <tr1/functional>
#include <tr1/memory>

#include <iostream>

using namespace honei;
using namespace tests;

extern "C"
{
    extern spe_program_handle_t spu_kernel_test;
}

// Currently disabled, as we don't yet support kernel switching!
#if 0
class SPEKernelTest :
    public QuickTest
{
    public:
        SPEKernelTest() :
            QuickTest("spe_kernel_test")
        {
        }

        virtual void run() const
        {
            Environment e __attribute__((aligned(16)));
            e.begin = 0x4000;
            e.end = 0x4000;

            std::tr1::shared_ptr<SPEKernel> kernel(new SPEKernel(spu_kernel_test, &e));

            unsigned long long retval_finished(0x0FEDCBA987654321ll);
            unsigned retval_dword(0x0);
            unsigned long long retval_qword(0x0ll);
            Instruction inst;

            inst.opcode = oc_test_instruction_finished;
            inst.c = &retval_finished;
            unsigned finished(kernel->enqueue(inst));

            inst.opcode = oc_test_result_dword;
            inst.c = &retval_dword;
            unsigned result_dword(kernel->enqueue(inst));

            inst.opcode = oc_test_result_qword;
            inst.c = &retval_qword;
            unsigned result_qword(kernel->enqueue(inst));

            inst.opcode = oc_halt;
            inst.c = 0;
            unsigned halt(kernel->enqueue(inst));

            std::cout << "XXX: before running!" << std::endl;
            SPEManager::instance()->begin()->run(*kernel);

            std::cout << "XXX: before first wait!" << std::endl;
            kernel->wait(finished, inst);
            TEST_CHECK_EQUAL(inst.opcode, oc_test_instruction_finished);
            TEST_CHECK_EQUAL(retval_finished, 0x0FEDCBA987654321);

            kernel->wait(result_dword, inst);
            TEST_CHECK_EQUAL(inst.opcode, oc_test_result_dword);
            TEST_CHECK_EQUAL(retval_dword, 0x12345678);

            std::cout << "XXX: waiting for qword!" << std::endl;
            kernel->wait(result_qword, inst);
            std::cout << "XXX: stopped waiting for qword!" << std::endl;
            TEST_CHECK_EQUAL(inst.opcode, oc_test_result_qword);
            TEST_CHECK_EQUAL(retval_qword, 0x123456789ABCDEF0);

            std::cout << "XXX: waiting for halt!" << std::endl;
            kernel->wait(halt, inst);
            std::cout << "XXX: halt arrived!" << std::endl;
        }
} spe_kernel_test;
#endif
