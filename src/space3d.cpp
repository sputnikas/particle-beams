#include "space3d.h"

//////////////////////////////////////////////////////////////////////////////
// Space3d implementation
//////////////////////////////////////////////////////////////////////////////

Space3d::~Space3d() {};

//////////////////////////////////////////////////////////////////////////////
// Space3dHalf implementation
//////////////////////////////////////////////////////////////////////////////

Space3dHalf::Space3dHalf() : normal(Vec3<double>(0,0,1)), point(Vec3<double>()) {}
Space3dHalf::Space3dHalf(Vec3<double> normal) : normal(normal/normal.norm()), point(Vec3<double>()) {}
Space3dHalf::Space3dHalf(Vec3<double> normal, Vec3<double> point) : normal(normal/normal.norm()), point(point) {}

int Space3dHalf::include(const Vec3<double> &r) {
    return (r - point).dot(normal) > 0;
}

//////////////////////////////////////////////////////////////////////////////
// Общие функции
//////////////////////////////////////////////////////////////////////////////

int in_spaces(const std::vector<Space3d*> &spaces, const Vec3<double> r) {
    for (std::vector<Space3d*>::const_iterator i = spaces.cbegin(); i != spaces.cend(); i++) {
        if ((*i)->include(r)) {
            return 1;
        }
    }
    return 0;
}