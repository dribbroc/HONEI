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
        enum CommModes
        {
            com_unidir_left_send_right_receive = 0,
            com_unidir_left_receive_right_send,
            com_bidir_send_receive,
            com_bidir_average
                //TODO...
        };

        template<unsigned i_ = 1, typename AttributeType_ = double, typename Tag_ = tags::CPU>
        class Communication
        {
            public:
                //example: Halo-based
                template<typename HaloType_>
                static void execute(HaloType_ & halo, unsigned attr_index, CommModes mode = com_bidir_send_receive) //TODO add rank
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
                static void execute(HaloType_ & halo, unsigned attr_index, CommModes mode = com_bidir_send_receive) //TODO add rank
                {
#ifdef FEM_COMM_DEBUG
                    if(halo.get_overlap() != 0)
                        throw CommunicationHaloOverlapMismatch(0, halo.get_overlap());
#endif
                }
                //TODO
        };

        template<typename Tag_>
        class Comm
        {
        };

        template<>
        class Comm<tags::CPU>
        {
            public:
                //example shared-mem exchange
                template<typename DataType_>
                static void send_receive(DataType_ * sendbuf, DataType_ * recvbuf, unsigned long num_elements, unsigned sender_rank = 0, unsigned recvr_rank = 0)
                {
                    DataType_ buf;
                    for(unsigned long i(0) ; i < num_elements ; ++i)
                    {
                        buf = sendbuf[i];
                        sendbuf[i] = recvbuf[i];
                        recvbuf[i] = sendbuf[i];
                    }
                }

                //TODO
        };
    }
}
#endif
