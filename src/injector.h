#pragma once

#include <vector>

#include "vec3.h"
#include "particle2.h"
#include "pmath.h"

//////////////////////////////////////////////////////////////////////////////
// Injector
// хранит данные по инжекторам
// type - номер типа, для данного типа вполне конкретен
//////////////////////////////////////////////////////////////////////////////

class Injector {
public:
    size_t type;
    
    Injector(const size_t &type);

    virtual ~Injector();

    virtual void inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt) = 0;
};

//////////////////////////////////////////////////////////////////////////////
// InjectorRectangleZ
// type  = 0
// Rectangular injector with homogenous in v and on S with Z direction of v
//////////////////////////////////////////////////////////////////////////////
class InjectorRectangleZ : public Injector {
public:
    Vec3<double> center;
    double wx;
    double wy;
    double vz;
    size_t N;
    ParticleType *ptype;
    double I;
    double rho;
    double kgeq;

    InjectorRectangleZ();
    InjectorRectangleZ( const double &wx, 
                        const double &wy, 
                        const double &vz, 
                        const size_t &N, 
                        ParticleType *ptype, 
                        const double &I, 
                        const double &rho, 
                        const double &kgeq);
    InjectorRectangleZ( const Vec3<double> &center,
                        const double &wx, 
                        const double &wy, 
                        const double &vz, 
                        const size_t &N, 
                        ParticleType *ptype, 
                        const double &I, 
                        const double &rho, 
                        const double &kgeq);

    void inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt);
};

//////////////////////////////////////////////////////////////////////////////
// InjectorRingZ
// type  = 1
// Ring injector with homogenous in v and on S with Z direction of v
//////////////////////////////////////////////////////////////////////////////
class InjectorRingZ : public Injector {
public:
    Vec3<double> center;
    double ra;
    double rb;
    double alpha0;
    double alpha;
    double vz;
    size_t N;
    ParticleType *ptype;
    double I;
    double rho;
    double kgeq;

    InjectorRingZ();
    InjectorRingZ(  const double &ra,
                    const double &rb,
                    const double &alpha0,
                    const double &alpha,
                    const double &vz, 
                    const size_t &N, 
                    ParticleType *ptype, 
                    const double &I, 
                    const double &rho, 
                    const double &kgeq);
    InjectorRingZ(  const Vec3<double> &center,
                    const double &ra,
                    const double &rb,
                    const double &alpha0,
                    const double &alpha,
                    const double &vz, 
                    const size_t &N, 
                    ParticleType *ptype, 
                    const double &I, 
                    const double &rho, 
                    const double &kgeq);

    void inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt);
};