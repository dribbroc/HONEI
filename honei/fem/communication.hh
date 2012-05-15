/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_COMMUNICATION_HH
#define FEM_GUARD_COMMUNICATION_HH 1

#include<honei/fem/halo.hh>
#include<honei/fem/communication_error.hh>
#include<vector>


namespace honei
{
    namespace fem
    {
        template<unsigned i_ = 1, typename AttributeType_ = double, typename Tag_ = tags::CPU>
        class Communication
        {
            public:
                template<typename HaloType_>
                static void execute(HaloType_ & halo, unsigned attr_index)
                {
#ifdef FEM_COMM_DEBUG
                    if(i_ != halo.get_overlap())
                        throw CommunicationHaloOverlapMismatch(i_, halo.get_overlap());
#endif
                }
                //TODO
        };

        template<typename AttributeType_, typename Tag_>
        class Communication<0, AttributeType_, Tag_>
        {
            public:
                template<typename HaloType_>
                static void execute(HaloType_ & halo, unsigned attr_index)
                {
#ifdef FEM_COMM_DEBUG
                    if(halo.get_overlap() != 0)
                        throw CommunicationHaloOverlapMismatch(0, halo.get_overlap());
#endif
                }
                //TODO
        };
    }
}
#endif
