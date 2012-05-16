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
            com_send = 0,
            com_receive,
            com_send_receive,
            com_average
                //TODO...
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
                static inline void send_recv(DataType_ * sendbuf, unsigned long dest_rank, unsigned long num_elements, DataType_* recvbuf, unsigned long source_rank)
                {
                    DataType_ buf;
                    const unsigned long i_end(num_elements);
                    for(unsigned long i(0) ; i < i_end ; ++i)
                    {
                        buf = sendbuf[i];
                        sendbuf[i] = recvbuf[i];
                        recvbuf[i] = buf;
                    }
                }

                //TODO
        };

        //overlap is a compile-time decision now, if not feasible, move to inner function template
        template<unsigned i_ = 1, CommModes j_ = com_send_receive, typename AttributeType_ = double, typename Tag_ = tags::CPU>
        class Communication
        {
            public:
                //example: Halo-based
                template<typename HaloType_, typename MeshType_>
                static void execute(HaloType_ & halo,
                                    unsigned attr_index,
                                    MeshType_ & other_mesh,
                                    unsigned long other_rank) //TODO other_mesh resolved by process mesh list (mesh by id), needs to be a local thing
                {
#ifdef FEM_COMM_DEBUG
                    if(i_ != halo.get_overlap())
                        throw CommunicationHaloOverlapMismatch(i_, halo.get_overlap());
#endif
                    //temp example
                    switch(j_)
                    {
                        case com_send_receive:
                            {
                                //temp until some sort of patch-internal buffer or master bufferpool available. Emulates Pack(). if rank=.. states here
                                AttributeType_* sendbuf(new AttributeType_[halo.size()]);
                                AttributeType_* recvbuf(new AttributeType_[halo.size()]);
                                for(unsigned long i(0) ; i < halo.size() ; ++i)
                                {
                                    if(typeid(AttributeType_) == typeid(typename MeshType_::attr_type_1_))
                                    {
                                        sendbuf[i] = halo.get_mesh().get_attributes_of_type_1().at(attr_index).at(halo.get_element(i)); //TODO maybe mixed up send/recv
                                        recvbuf[i] = other_mesh.get_attributes_of_type_1().at(attr_index).at(halo.get_element_counterpart(i));
                                    }
                                }
                                //'post'
                                fem::Comm<Tag_>::send_recv(sendbuf, other_rank, halo.size(), recvbuf, halo.get_mesh().get_pp_rank());
                                for(unsigned long i(0) ; i < halo.size() ; ++i)
                                {
                                    if(typeid(AttributeType_) == typeid(typename MeshType_::attr_type_1_))
                                    {
                                        halo.get_mesh().get_attributes_of_type_1().at(attr_index).at(halo.get_element(i)) = sendbuf[i]; //TODO maybe mixed up send/recv
                                        other_mesh.get_attributes_of_type_1().at(attr_index).at(halo.get_element_counterpart(i)) = recvbuf[i];
                                    }
                                }
                            }
                    }
                }
                //TODO
        };

    }
}
#endif
