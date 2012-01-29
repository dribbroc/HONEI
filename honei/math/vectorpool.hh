/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MATH_GUARD_VECTORPOOL_HH
#define MATH_GUARD_VECTORPOOL_HH 1

#include<vector>

namespace honei
{
    template<typename VT_>
    static std::vector<VT_>  create_vectorpool(unsigned long pool_size, unsigned long vector_size)
    {
        std::vector<VT_> result;
        for(unsigned long i(0) ; i < pool_size ; ++i)
        {
            VT_ vector(vector_size);
            result.push_back(vector);
        }
        return result;
    }
}
#endif
