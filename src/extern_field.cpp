#include "extern_field.h"

//////////////////////////////////////////////////////////////////////////////
// ExternField - abstract
//////////////////////////////////////////////////////////////////////////////

ExternField::~ExternField() {};

//////////////////////////////////////////////////////////////////////////////
// ExternFieldConst - постоянные поля
//////////////////////////////////////////////////////////////////////////////

ExternFieldConst::ExternFieldConst(const Vec3<double> &E, const Vec3<double> &B) : E(E), B(B) {}
ExternFieldConst::ExternFieldConst(const Vec3<double> &E, const Vec3<double> &B, Space3d *space3d) : E(E), B(B) { spaces.push_back(space3d); }
    
void ExternFieldConst::field(Vec3<double> &Ee, Vec3<double> &Be, const Vec3<double> &r, const double &t) {
    if (in_spaces(spaces, r)) {
        Ee += E;
        Be += B;
    }
}