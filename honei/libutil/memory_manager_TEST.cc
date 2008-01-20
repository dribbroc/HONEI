/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/assertion.hh>
#include <honei/libutil/memory_manager.hh>
#include <honei/libutil/memory_backend.hh>
#include <honei/libutil/stringify.hh>
#include <unittest/unittest.hh>

#include <cmath>
#include <map>
#include <iostream>

using namespace honei;
using namespace tests;

namespace
{
    class FakeBackend :
        public MemoryBackend
    {
        private:
            struct Info
            {
                char * address;
                unsigned long size;

                Info(char * a, unsigned long s) :
                    address(a),
                    size(s)
                {
                }
            };

            std::map<MemoryId, Info *> _info_map;

        public:
            static FakeBackend * instance()
            {
                static FakeBackend result;

                return &result;
            }

            static MemoryBackend * backend_instance()
            {
                return instance();
            }

            virtual void upload(const MemoryId id, const DeviceId, void * address, const std::ptrdiff_t size)
            {
                const char * a(static_cast<const char *>(address));
                Info * info(new Info(new char[size], size));

                std::copy(a, a + size, info->address);

                _info_map.insert(std::make_pair(id, info));
            }

            virtual void download(const MemoryId id, const DeviceId, void * address, const std::ptrdiff_t size)
            {
                char * a(static_cast<char *>(address));
                std::map<MemoryId, Info *>::iterator i(_info_map.find(id));
                ASSERT(_info_map.end() != i, "Unknown id '" + stringify(id) + "'!");
                ASSERT(size <= i->second->size, "Tried to download more memory than was uploaded!");

                std::copy(i->second->address, i->second->address + size, a);
            }

            virtual void free(const MemoryId id, const DeviceId)
            {
                std::map<MemoryId, Info *>::iterator i(_info_map.find(id));
                ASSERT(_info_map.end() != i, "Unknown id '" + stringify(id) + "'!");

                delete[] i->second->address;

                _info_map.erase(i);
            }

            virtual void swap(const MemoryId left, const MemoryId right)
            {
                std::map<MemoryId, Info *>::iterator l(_info_map.find(left)),
                    r(_info_map.find(right));
                ASSERT(_info_map.end() != l, "Unknown id '" + stringify(left) + "'!");
                ASSERT(_info_map.end() != r, "Unknown id '" + stringify(right) + "'!");

                Info * tmp(l->second);
                l->second = r->second;
                r->second = tmp;
            }
    };
}

static MemoryBackendRegistrator fake_backend_registrator(tags::tv_fake, &FakeBackend::backend_instance);

class MemoryManagerTest :
    public BaseTest
{
    public:
        MemoryManagerTest() :
            BaseTest("memory_manager_test")
        {
        }

        virtual void run() const
        {
            float a[30], b[30];

            std::fill_n(a, 30, float(1.11));
            std::fill_n(b, 30, float(2.22));

            MemoryId id_a(MemoryManager::instance()->associate(a, 30 * sizeof(float))),
                    id_b(MemoryManager::instance()->associate(b, 30 * sizeof(float)));
            MemoryManager::instance()->upload(id_a, default_device, tags::tv_fake, 30 * sizeof(float));
            MemoryManager::instance()->upload(id_b, default_device, tags::tv_fake, 30 * sizeof(float));

            MemoryManager::instance()->download(id_b, default_device);
            bool ok(true);
            for (float * i(b) ; i < b + 30 ; ++i)
            {
                ok &= (fabs(*i - 2.22) <= std::numeric_limits<float>::epsilon());
            }
            TEST_CHECK(ok);

            MemoryManager::instance()->download(id_a, default_device);
            ok = true;
            for (float * i(a) ; i < a + 30 ; ++i)
            {
                ok &= (fabs(*i - 1.11) <= std::numeric_limits<float>::epsilon());
            }
            TEST_CHECK(ok);

            MemoryManager::instance()->swap(id_a, id_b);

            MemoryManager::instance()->download(id_b, default_device);
            ok = true;
            for (float * i(a) ; i < a + 30 ; ++i)
            {
                ok &= (fabs(*i - 1.11) <= std::numeric_limits<float>::epsilon());
            }
            TEST_CHECK(ok);

            MemoryManager::instance()->download(id_a, default_device);
            ok = true;
            for (float * i(b) ; i < b + 30 ; ++i)
            {
                ok &= (fabs(*i - 2.22) <= std::numeric_limits<float>::epsilon());
            }
            TEST_CHECK(ok);
        }
} memory_manager_test;
