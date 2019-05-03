#include "ppsolver.h"

int PPSolver::in_spaces(const Vec3<double> r1, const Vec3<double> r2) {
    for (std::vector<Space3d*>::iterator i = spaces.begin(); i != spaces.end(); i++) {
        if ((*i)->include(r1) && (*i)->include(r2)) {
            return 1;
        }
    }
    return 0;
}

int PPSolver::field_interaction1(Vec3<double> &Ei, Vec3<double> &Bi, const std::vector<Particle*>::iterator &i) {
    Ei = Vec3<double>();
    Bi = Vec3<double>();
    for (std::vector<Particle*>::iterator j = particles.begin(); j != particles.end(); j++) {
        if (j != i) {
            std::pair<size_t, ParticleData> result;
            Vec3<double> r;
            std::vector<ParticleData>::reverse_iterator pj;
            switch (interaction_type) {
            case 0:
                switch (sphere_method) {
                case 0:
                    r = (*i)->p.rbegin()->r;
                    result = sm_bisection(*j, r, t, dt);
                    if ((result.first != (size_t) -1) && in_spaces(r, result.second.r)) {
                        (*i)->type->field(Ei, Bi, result.second, r);
                    }
                    break;
                case 1:
                    r = (*i)->p.rbegin()->r;
                    result = sm_left_bisection(*j, r, t, dt);
                    if ((result.first != (size_t) -1) && in_spaces(r, result.second.r)) {
                        (*i)->type->field(Ei, Bi, result.second, r);
                    }
                    break;
                case 2:
                    r = (*i)->p.rbegin()->r;
                    result = sm_left_intervals(*j, r, t, dt);
                    if ((result.first != (size_t) -1) && in_spaces(r, result.second.r)) {
                        (*i)->type->field(Ei, Bi, result.second, r);
                    }
                    break;
                case 3:
                    r = (*i)->p.rbegin()->r;
                    pj = (*j)->p.rbegin();
                    if (in_spaces(r, pj->r)) {
                        (*i)->type->field(Ei, Bi, *pj, r);
                    }
                    break;
                }
                break;
            case 1:
                switch (sphere_method) {
                case 0:
                    r = (*i)->p.rbegin()->r;
                    pj = (*j)->p.rbegin();
                    if (in_spaces(r, pj->r)) {
                        (*i)->type->fieldN(Ei, Bi, *pj, r);
                    }
                    break;
                }
                break;
            }
        }
    }
    return 0;
}

int PPSolver::field_extern(Vec3<double> &Ee, Vec3<double> &Be, const Vec3<double> &r, const double &t) {
    Ee = Vec3<double>();
    Be = Vec3<double>();
    for (std::vector<ExternField*>::iterator j = extern_fields.begin(); j != extern_fields.end(); j++) {
        (*j)->field(Ee, Be, r, t);
    }
    return 0;
}

ParticleData PPSolver::euler(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    field_interaction(Ei, Bi, i);
    field_extern(Ee, Be, pi->r, t);
    pi->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v), pi->v);
    p.v = pi->v + pi->a*CL*dt;
    p.r = pi->r + pi->v*CL*dt;
    return p;
}

ParticleData PPSolver::euler_u(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    field_interaction(Ei, Bi, i);
    field_extern(Ee, Be, pi->r, t);
    Vec3<double> fpm = (*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v);
    pi->a = acceleration(fpm, pi->v);
    p.v = velocity_u(momentum_v(pi->v) + fpm*CL*dt);
    p.r = pi->r + pi->v*CL*dt;
    return p;
}

ParticleData PPSolver::runge4v_mistake(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be, k1r, k2r, k3r, k4r, k1v, k2v, k3v, k4v;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    field_interaction(Ei, Bi, i);
    
    field_extern(Ee, Be, pi->r, t);
    pi->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v), pi->v);
    k1r = pi->v*cdt;
    k1v = pi->a*cdt;

    field_extern(Ee, Be, pi->r + 0.5*k1r, t + dt/2);
    k2r = (pi->v + 0.5*k1v)*cdt;
    k2v = (acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v + 0.5*k1v), pi->v + 0.5*k1v))*cdt;

    field_extern(Ee, Be, pi->r + 0.5*k2r, t + dt/2);
    k3r = (pi->v + 0.5*k2v)*cdt;
    k3v = (acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v + 0.5*k2v), pi->v + 0.5*k2v))*cdt;

    field_extern(Ee, Be, pi->r + k3r, t + dt);
    k4r = (pi->v + k3v)*cdt;
    k4v = (acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v + k3v), pi->v + k3v))*cdt;

    p.v = p.v + 1./6*(k1v + 2.*k2v + 2.*k3v + k4v);
    p.r = p.r + 1./6*(k1r + 2.*k2r + 2.*k3r + k4r);
    return p;
}

ParticleData PPSolver::runge4u_mistake(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be, k1r, k2r, k3r, k4r, k1u, k2u, k3u, k4u;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    Vec3<double> pu = momentum_v(pi->v);
    field_interaction(Ei, Bi, i);
    
    field_extern(Ee, Be, pi->r, t);
    pi->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v), pi->v);
    k1r = pi->v*cdt;
    k1u = (*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v)*cdt;

    field_extern(Ee, Be, pi->r + 0.5*k1r, t + dt/2);
    k2r = velocity_u(pu + 0.5*k1u)*cdt;
    k2u = (*i)->type->forcePm(Ee + Ei, Be + Bi, velocity_u(pu + 0.5*k1u))*cdt;

    field_extern(Ee, Be, pi->r + 0.5*k2r, t + dt/2);
    k3r = velocity_u(pu + 0.5*k2u)*cdt;;
    k3u = (*i)->type->forcePm(Ee + Ei, Be + Bi, velocity_u(pu + 0.5*k2u))*cdt;

    field_extern(Ee, Be, pi->r + k3r, t + dt);
    k4r = velocity_u(pu + k3u)*cdt;;
    k4u = (*i)->type->forcePm(Ee + Ei, Be + Bi, velocity_u(pu + k3u))*cdt;

    p.v = velocity_u(pu + 1./6*(k1u + 2.*k2u + 2.*k3u + k4u));
    p.r = p.r + 1./6*(k1r + 2.*k2r + 2.*k3r + k4r);
    return p;
}

ParticleData PPSolver::adams_bashfort4_v(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator p1 = (*i)->p.rbegin();
    std::vector<ParticleData>::reverse_iterator p2 = p1++;
    std::vector<ParticleData>::reverse_iterator p3 = p2++;
    std::vector<ParticleData>::reverse_iterator p4 = p3++;
    std::vector<ParticleData>::reverse_iterator rend = (*i)->p.rend();
    ParticleData p;

    field_interaction(Ei, Bi, i);
    field_extern(Ee, Be, p1->r, t);
    p1->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, p1->v), p1->v);
    
    if (p4 != rend) {
        p.v = p1->v + cdt*(55./24.*p1->a - 59./24.*p2->a + 37./24.*p3->a - 9./24.*p4->a);
        p.r = p1->r + cdt*(55./24.*p1->v - 59./24.*p2->v + 37./24.*p3->v - 9./24.*p4->v);
    } else if (p3 != rend) {
        p.v = p1->v + cdt*(23./12.*p1->a - 16./12.*p2->a + 5./12.*p3->a);
        p.r = p1->r + cdt*(23./12.*p1->v - 16./12.*p2->v + 5./12.*p3->v);
    } else if (p2 != rend) {
        p.v = p1->v + cdt*(1.5*p1->a - 0.5*p2->a);
        p.r = p1->r + cdt*(1.5*p1->v - 0.5*p2->v);
    } else {
        p.v = p1->v + cdt*p1->a;
        p.r = p1->r + cdt*p1->v;
    }
    return p;
}

ParticleData PPSolver::adams_bashfort4_u(const std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator p1 = (*i)->p.rbegin();
    std::vector<ParticleData>::reverse_iterator p2 = p1++;
    std::vector<ParticleData>::reverse_iterator p3 = p2++;
    std::vector<ParticleData>::reverse_iterator p4 = p3++;
    std::vector<ParticleData>::reverse_iterator rend = (*i)->p.rend();
    ParticleData p;

    field_interaction(Ei, Bi, i);
    field_extern(Ee, Be, p1->r, t);
    p1->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, p1->v), p1->v);
    
    if (p4 != rend) {
        p.v = velocity_u(momentum_v(p1->v) 
            + cdt*(55./24.*forcepm(p1->a, p1->v) 
            - 59./24.*forcepm(p2->a, p2->v) 
            + 37./24.*forcepm(p3->a, p3->v) 
            -  9./24.*forcepm(p4->a, p4->v)));
        p.r = p1->r + cdt*(55./24.*p1->v - 59./24.*p2->v + 37./24.*p3->v - 9./24.*p4->v);
    } else if (p3 != rend) {
        p.v = velocity_u(momentum_v(p1->v) 
            + cdt*(23./12.*forcepm(p1->a, p1->v) 
            - 16./12.*forcepm(p2->a, p2->v) 
            +  5./12.*forcepm(p3->a, p3->v)));
        p.r = p1->r + cdt*(23./12.*p1->v - 16./12.*p2->v + 5./12.*p3->v);
    } else if (p2 != rend) {
        p.v = velocity_u(momentum_v(p1->v) 
            + cdt*(1.5*forcepm(p1->a, p1->v) 
            - 0.5*forcepm(p2->a, p2->v)));
        p.r = p1->r + cdt*(1.5*p1->v - 0.5*p2->v);
    } else {
        p.v = velocity_u(momentum_v(p1->v) + cdt*forcepm(p1->a, p1->v));
        p.r = p1->r + cdt*p1->v;
    }
    return p;
}

///////////////////////////////////////////////////////////////////////////////
// Основной расчётный цикл
///////////////////////////////////////////////////////////////////////////////

size_t PPSolver::calc() {
    nt++;
    t += dt;
    ct += cdt;
    for (std::vector<Injector*>::iterator i = injectors.begin(); i != injectors.end(); i++) {
        (*i)->inject(particles, cdt, nt);
    }
    if (number_method < RUNGE4V) {
        for (std::vector<Particle*>::iterator i = particles.begin(); i != particles.end() && (*i)->nmax == nt - 1; i++) {
            ParticleData p;
            switch (number_method) {
            case EULER:
                p = euler(i);
                break;
            case EULER_U:
                p = euler_u(i);
                break;
            case RUNGE4V_MISTAKE:
                p = runge4v_mistake(i);
                break;
            case RUNGE4U_MISTAKE:
                p = runge4u_mistake(i);
                break;
            case ADAMS_BASHFORT4_V:
                p = adams_bashfort4_v(i);
                break;
            case ADAMS_BASHFORT4_U:
                p = adams_bashfort4_u(i);
                break;
            }
            int inspace = 0;
            for (std::vector<Space3d*>::iterator j = spaces.begin(); j != spaces.end(); j++) {
                inspace += (*j)->include(p.r);
            }
            if (inspace != 0) {
                (*i)->nmax = nt;
                (*i)->p.push_back(p);
            }
        }
    } else {

    }
    return nt;
}