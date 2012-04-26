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
            typename StorageType_ = std::vector<IndexType_>,
            template<typename, typename> class OuterAttrStorageType_ = std::vector,
            template<typename, typename> class AttrStorageType_ = std::vector,
            typename AttributeType1_ = Nil,
            typename AttributeType2_ = Nil,
            typename AttributeType3_ = Nil,
            typename AttributeType4_ = Nil>
        class Topology
        {
            public:
                Topology() :
                    _num_polytopes(0),
                    _num_attributes(0),
                    _topology(new OuterStorageType_<StorageType_, std::allocator<StorageType_> >),
                    _attributes_of_type_1(new OuterAttrStorageType_< //todo just do it if needed
                    AttrStorageType_<AttributeType1_, std::allocator<AttributeType1_> >,
                    std::allocator<AttrStorageType_<AttributeType1_, std::allocator<AttributeType1_> > > >)
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
                unsigned long _num_attributes;
                OuterStorageType_<StorageType_, std::allocator<StorageType_> >* _topology;
                OuterAttrStorageType_<
                    AttrStorageType_<AttributeType1_, std::allocator<AttributeType1_> >, //todo for other types
                    std::allocator<AttrStorageType_<AttributeType1_, std::allocator<AttributeType1_> > > >* _attributes_of_type_1;
        };
    }
}
#endif
