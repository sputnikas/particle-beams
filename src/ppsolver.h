#pragma once

#include <vector>

#include "vec3.h"
#include "particle2.h"
#include "injector.h"
#include "extern_field.h"
#include "space3d.h"

class PPSolver {
public:
    enum NUMBER_METHOD {
        EULER,
        EULER_U,
        RUNGE4V_MISTAKE,
        RUNGE4U_MISTAKE,
        ADAMS_BASHFORT4_V,
        ADAMS_BASHFORT4_U,
        RUNGE4V
    };
    std::vector<Space3d*>      spaces;
    std::vector<Particle*>     particles;
    std::vector<ParticleType*> particle_types;
    std::vector<Injector*>     injectors;
    std::vector<ExternField*>  extern_fields;
    size_t (*field_interaction)(Vec3<double> &Ei, Vec3<double> &Bi, const std::vector<Particle*>::iterator &i);
    size_t interaction_type;
    size_t sphere_method;
    size_t number_method;
    size_t nt;
    double dt;
    double t;
    double cdt;
    double ct;

    int in_spaces(const Vec3<double> r1, const Vec3<double> r2);
    int field_interaction1(Vec3<double> &Ei, Vec3<double> &Bi, const std::vector<Particle*>::iterator &i);
    int field_extern(Vec3<double> &Ee, Vec3<double> &Be, const Vec3<double> &r, const double &t);
    ParticleData euler(const std::vector<Particle*>::iterator &i);
    ParticleData euler_u(const std::vector<Particle*>::iterator &i);
    ParticleData runge4v_mistake(const std::vector<Particle*>::iterator &i);
    ParticleData runge4u_mistake(const std::vector<Particle*>::iterator &i);
    ParticleData adams_bashfort4_v(const std::vector<Particle*>::iterator &i);
    ParticleData adams_bashfort4_u(const std::vector<Particle*>::iterator &i);
    size_t calc();
};