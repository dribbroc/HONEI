/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_TOPOLOGY_HH
#define FEM_GUARD_TOPOLOGY_HH 1

#include <vector>

namespace honei
{
    class Nil;

    namespace fem
    {
        //For compatibility: storage types must only have size(), push_back() and [] members in common
        template<
            typename IndexType_ = unsigned long,
            template<typename, typename> class OuterStorageType_ = std::vector,
            typename StorageType_ = std::vector<IndexType_> >
        class Topology
        {
            public:
                typedef IndexType_ index_type_;
                typedef StorageType_ storage_type_;
                typedef OuterStorageType_<StorageType_, std::allocator<StorageType_> > outer_storage_type_;

                Topology() :
                    _num_polytopes(0),
                    _topology(new OuterStorageType_<StorageType_, std::allocator<StorageType_> >)
                {
                };

                ~Topology()
                {
                    delete _topology;
                }

                unsigned long size()
                {
                    return _num_polytopes;
                }

                void push_back(const StorageType_ s)
                {
                    _topology->push_back(s);
                    ++_num_polytopes;
                }

                StorageType_ & at(unsigned long i)
                {
                    return _topology->at(i);
                }

                StorageType_ & operator[] (unsigned long i)
                {
                    return _topology->at(i);
                }

                void push_back()
                {
                    StorageType_ s;
                    _topology->push_back(s);
                    ++_num_polytopes;
                }

            private:
                unsigned long _num_polytopes;
                OuterStorageType_<StorageType_, std::allocator<StorageType_> >* _topology;
        };

    }
}
#endif
