#include <iostream>
#include <vector>

#include "vec3.h"

const double EQ = 1.6021892E-19;
const double EM = 9.109534E-31;
const double EPS0 = 8.854187817E-12;
const double K0 = 8.987551788E9;
const double CL = 2.99792458E8;
const double RL = 2.8179402894E-15;
const double CL2 = 8.9875517873681764E16;

//////////////////////////////////////////////////////////////////////////////
// Particle
// После долгих размышлений появилась одна идея:
//      частицы, должны быть упорядочены во временной массив!
// 
//////////////////////////////////////////////////////////////////////////////

class ParticleTypeBase {
public:
    
};

template <typename PType>
class ParticleDataBase {

};

class ParticleTypePointR : public ParticleTypeBase {
public:
    double q;
    double m;
    double qPm;
};

class ParticleTypeSphereR : public ParticleTypeBase {
public:
    double q;
    double m;
    double qPm;
    double radius;
};

class ParticleDataBase<ParticleTypePointR> {
public:
    Vec3<double> r;
    Vec3<double> v;
    Vec3<double> a;

    void field(Vec3<double> &E, Vec3<double> &B, const Vec3<double> &r, ParticleTypePointR *type);
};

void ParticleDataBase<ParticleTypePointR>::field(Vec3<double> &E, Vec3<double> &B, const Vec3<double> &rs, ParticleTypePointR *type) {
    Vec3<double> Rv = rs - r;
    double R = norm(Rv);
    Vec3<double> normal = Rv/R;

}

class ParticleBase {
public:

};

template <typename PType>
class Particle : public ParticleBase {
public:
    std::vector<ParticleDataBase<PType>> p;
    PType *type;
};