#pragma once

#include <vector>

#include "vec3.h"
#include "space3d.h"

//////////////////////////////////////////////////////////////////////////////
// ExternField - abstract
//////////////////////////////////////////////////////////////////////////////

class ExternField {
public:
    virtual ~ExternField();

    virtual void field(Vec3<double> &Ee, Vec3<double> &Be, const Vec3<double> &r, const double &t) = 0;
};

//////////////////////////////////////////////////////////////////////////////
// ExternFieldConst - постоянные однородные поля
//////////////////////////////////////////////////////////////////////////////

class ExternFieldConst : public ExternField {
public:
    std::vector<Space3d*> spaces;
    Vec3<double> E;
    Vec3<double> B;

    ExternFieldConst(const Vec3<double> &E, const Vec3<double> &B);
    ExternFieldConst(const Vec3<double> &E, const Vec3<double> &B, Space3d *space3d);
    void field(Vec3<double> &Ee, Vec3<double> &Be, const Vec3<double> &r, const double &t);
};