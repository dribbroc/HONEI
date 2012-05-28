/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef QT_GUARD_NORMAL_HH
#define QT_GUARD_NORMAL_HH 1

#include<cmath>
#include<clients/qt/vertex.hh>


template<typename DT_>
void get_unit_normal(Vertex<DT_> & v0, Vertex<DT_> & v1, Vertex<DT_> & v2, Vertex<DT_> & normal)
{
    Vertex<DT_> a, b;
    a.coord_x = v0.coord_x - v1.coord_x;
    a.coord_y = v0.coord_y - v1.coord_y;
    a.coord_z = v0.coord_z - v1.coord_z;

    b.coord_x = v1.coord_x - v2.coord_x;
    b.coord_y = v1.coord_y - v2.coord_y;
    b.coord_z = v1.coord_z - v2.coord_z;

    normal.coord_x = (a.coord_y * b.coord_z) - (a.coord_z * b.coord_y);
    normal.coord_y = (a.coord_z * b.coord_x) - (a.coord_x * b.coord_z);
    normal.coord_z = (a.coord_x * b.coord_y) - (a.coord_y * b.coord_x);

    DT_ len = (DT_)(sqrt((normal.coord_x * normal.coord_x) + (normal.coord_y * normal.coord_y) + (normal.coord_z * normal.coord_z)));

    len = len < std::numeric_limits<DT_>::epsilon() ? DT_(1) : len;

    normal.coord_x /= len;
    normal.coord_y /= len;
    normal.coord_z /= len;
}

#endif
