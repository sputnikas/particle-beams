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
    ParticleData p;
    p.r.x = center.x - 0.5*wx + random()*wx;
    p.r.y = center.y - 0.5*wy + random()*wy;
    p.r.z = center.z + random()*vz*cdt;
    p.v.x = 0;
    p.v.y = 0;
    p.v.z = vz;
    p.a.x = 0;
    p.a.y = 0;
    p.a.z = 0;
    Particle* result;
    result->p.push_back(p);
    result->type = ptype;
    result->nmax = nt;
    particles.push_back(result);
};

//////////////////////////////////////////////////////////////////////////////
// InjectorRingZ Implementation
//////////////////////////////////////////////////////////////////////////////

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
    Injector(0), center(Vec3<double>()), ra(ra), rb(rb), alpha0(alpha0), alpha(alpha), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
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
    Injector(0), center(center), ra(ra), rb(rb), alpha0(alpha0), alpha(alpha), vz(vz), N(N), ptype(ptype), I(I), rho(rho), kgeq(kgeq)
{
}; 

void InjectorRingZ::inject(std::vector<Particle*> &particles, const double &cdt, const size_t &nt) {
    ParticleData p;
    double R = ra + (rb - ra)*random();
    double phi = alpha0 + (alpha - alpha0)*random();
    p.r.x = center.x + R*cos(phi);
    p.r.y = center.y + R*sin(phi);
    p.r.z = center.z + random()*vz*cdt;
    p.v.x = 0;
    p.v.y = 0;
    p.v.z = vz;
    p.a.x = 0;
    p.a.y = 0;
    p.a.z = 0;
    Particle* result;
    result->p.push_back(p);
    result->type = ptype;
    result->nmax = nt;
    particles.push_back(result);
};