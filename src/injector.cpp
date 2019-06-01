#include "injector.h"

//////////////////////////////////////////////////////////////////////////////
// Injector Implementation
//////////////////////////////////////////////////////////////////////////////

Injector::Injector(const size_t &type) : type(type) {};
Injector::~Injector() {};

//////////////////////////////////////////////////////////////////////////////
// InjectorRectangleZ Implementation
//////////////////////////////////////////////////////////////////////////////

InjectorRectangleZ::InjectorRectangleZ( const double &wx, 
                                        const double &wy, 
                                        const double &vz, 
                                        const size_t &N, 
                                        ParticleType *ptype, 
                                        const double &I, 
                                        const double &rho, 
                                        const double &kgeq) :
    Injector(0), center(Vec3<double>()), wx(wx), wy(wy), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
{
};
InjectorRectangleZ::InjectorRectangleZ( const Vec3<double> &center,
                                        const double &wx, 
                                        const double &wy, 
                                        const double &vz, 
                                        const size_t &N, 
                                        ParticleType *ptype, 
                                        const double &I, 
                                        const double &rho, 
                                        const double &kgeq) :
    Injector(0), center(center), wx(wx), wy(wy), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
{
}; 

void InjectorRectangleZ::inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt) {
    for (size_t i = 0; i<N; i++) {
        ParticleData p;
        p.r.x = center.x - 0.5*wx + pb::random()*wx;
        p.r.y = center.y - 0.5*wy + pb::random()*wy;
        p.r.z = center.z + pb::random()*vz*cdt;
        p.v.x = 0;
        p.v.y = 0;
        p.v.z = vz;
        p.a.x = 0;
        p.a.y = 0;
        p.a.z = 0;
        Particle* result = new Particle();
        //result->p = std::vector<ParticleData>();
        result->p.push_back(p);
        result->type = ptype;
        result->nmax = nt;
        particles.push_back(result);
    }
};

//////////////////////////////////////////////////////////////////////////////
// InjectorRingZ Implementation
//////////////////////////////////////////////////////////////////////////////

InjectorRingZ::InjectorRingZ() :
    Injector(1), center(Vec3<double>()), ra(0), rb(0), alpha0(0), alpha(0), vz(0), N(0), ptype(NULL), I(0), rho(0), kgeq(0)
{
};

InjectorRingZ::InjectorRingZ( const double &ra,
                              const double &rb,
                              const double &alpha0,
                              const double &alpha, 
                              const double &vz, 
                              const size_t &N, 
                              ParticleType *ptype, 
                              const double &I, 
                              const double &rho, 
                              const double &kgeq) :
    Injector(1), center(Vec3<double>()), ra(ra), rb(rb), alpha0(alpha0), alpha(alpha), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
{
};

InjectorRingZ::InjectorRingZ( const Vec3<double> &center,
                              const double &ra,
                              const double &rb,
                              const double &alpha0,
                              const double &alpha, 
                              const double &vz, 
                              const size_t &N, 
                              ParticleType *ptype, 
                              const double &I, 
                              const double &rho, 
                              const double &kgeq) :
    Injector(1), center(center), ra(ra), rb(rb), alpha0(alpha0), alpha(alpha), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
{
}; 

void InjectorRingZ::inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt) {
    for (size_t i = 0; i<N; i++) {
        ParticleData p;
        double R = ra + (rb - ra)*pb::random();
        double phi = alpha0 + (alpha - alpha0)*pb::random();
        p.r.x = center.x + R*cos(phi);
        p.r.y = center.y + R*sin(phi);
        p.r.z = center.z + pb::random()*vz*cdt;
        p.v.x = 0;
        p.v.y = 0;
        p.v.z = vz;
        p.a.x = 0;
        p.a.y = 0;
        p.a.z = 0;
        Particle* result = new Particle();
        result->p = std::vector<ParticleData>();
        result->p.push_back(p);
        result->type = ptype;
        result->nmax = nt;
        particles.push_back(result);
    }
};