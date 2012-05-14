/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_HALO_HH
#define FEM_GUARD_HALO_HH 1

#include<vector>

namespace honei
{
    namespace fem
    {
        template<unsigned delta = 1, template<typename, typename> class StorageType_ = std::vector, typename IndexType_ = unsigned long>
        class Halo
        {
            public:
                typedef IndexType_ index_type_;

                Halo() :
                    _halo_element_counterparts(new StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _halo_element_indices(new StorageType_<IndexType_, std::allocator<IndexType_> >())
                {
                }

                ~Halo()
                {
                    delete _halo_element_counterparts;
                    delete _halo_element_indices;
                }

                ///Add correspondence of i in mesh 1
                void add_halo_element_counterpart(IndexType_ i, IndexType_ j)
                {
                }

                IndexType_ get_halo_element_counterpart(IndexType_ i)
                {
                    return _halo_element_counterparts->at(_halo_element_indices->at(i));
                }

                template<typename MeshType_>
                IndexType_ get_element_counterpart(MeshType_ & source, MeshType_ & target, IndexType_ index)
                {
                    return 47;
                }

            private:
                StorageType_<IndexType_, std::allocator<IndexType_> >* _halo_element_counterparts;
                StorageType_<IndexType_, std::allocator<IndexType_> >* _halo_element_indices;
        };

        template<template<typename, typename> class StorageType_, typename IndexType_>
        class Halo<0u, StorageType_, IndexType_>
        {
            public:
                typedef IndexType_ index_type_;

                Halo() :
                    _halo_element_counterparts(new StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _halo_element_indices(new StorageType_<IndexType_, std::allocator<IndexType_> >())
                {
                }

                ~Halo()
                {
                    delete _halo_element_counterparts;
                    delete _halo_element_indices;
                }

                ///Add correspondence of i in mesh 1
                void add_halo_element_counterpart(IndexType_ i, IndexType_ j)
                {
                }

                IndexType_ get_halo_element_counterpart(IndexType_ i)
                {
                    return _halo_element_counterparts->at(_halo_element_indices->at(i));
                }

                template<typename MeshType_>
                IndexType_ get_element_counterpart(MeshType_ & source, MeshType_ & target, IndexType_ index)
                {
                    return 48;
                }

            private:
                StorageType_<IndexType_, std::allocator<IndexType_> >* _halo_element_counterparts;
                StorageType_<IndexType_, std::allocator<IndexType_> >* _halo_element_indices;
        };

        /*template<typename HaloType_, typename MeshType_, typename IndexType_>
        static IndexType_ get_element_counterpart(MeshType_ & source, MeshType_ & target, HaloType_ & halo, IndexType_ index)
        {
        }

        template<typename MeshType_, typename IndexType_>
        static IndexType_  get_element_counterpart(MeshType_ & source, MeshType_ & target, Halo<0, std::vector> & halo, IndexType_ index)
        {
            return 48;
        }*/

        /*template<typename IndexType_>
        static StorageType_ & get_primary_comm_neighbours(MeshType_ & m1, MeshType_ & m2, Halo<0, StorageType_, IndexType_> & halo, IndexType_ index)
        {
        }

        template<typename IndexType_>
        static StorageType_ & get_all_comm_neighbours(MeshType_ & m1, MeshType_ & m2, Halo<0, StorageType_, IndexType_> & halo, IndexType_ index)
        {
        }*/
    }

}


#endif
