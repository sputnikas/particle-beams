#include "ppsolver.h"

std::ofstream fd;

int PPSolver::in_spaces(const Vec3<double> r1, const Vec3<double> r2) {
    for (std::vector<Space3d*>::iterator i = spaces.begin(); i != spaces.end(); i++) {
        if ((*i)->include(r1) && (*i)->include(r2)) {
            return 1;
        }
    }
    return 0;
}

int PPSolver::field_interaction0(Vec3<double> &Ei, Vec3<double> &Bi, const std::vector<Particle*>::iterator &i) {
    Ei = Vec3<double>();
    Bi = Vec3<double>();
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

ParticleData PPSolver::euler(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    (this->*field_interaction)(Ei, Bi, i);
    field_extern(Ee, Be, pi->r, t);
    pi->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v), pi->v);
    p.v = pi->v + pi->a*CL*dt;
    p.r = pi->r + pi->v*CL*dt;
    //std::cout << pi->r << pi->v << pi->a << std::endl;
    //std::cout << Ei << Bi << std::endl;
    //getchar();
    return p;
}

ParticleData PPSolver::euler_u(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    (this->*field_interaction)(Ei, Bi, i);
    field_extern(Ee, Be, pi->r, t);
    Vec3<double> fpm = (*i)->type->forcePm(Ee + Ei, Be + Bi, pi->v);
    pi->a = acceleration(fpm, pi->v);
    p.v = velocity_u(momentum_v(pi->v) + fpm*CL*dt);
    p.r = pi->r + pi->v*CL*dt;
    return p;
}

ParticleData PPSolver::runge4v_mistake(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be, k1r, k2r, k3r, k4r, k1v, k2v, k3v, k4v;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    (this->*field_interaction)(Ei, Bi, i);
    
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

    //std::cout << k1r << k1v << std::endl;
    //std::cout << k2r << k2v << std::endl;
    //std::cout << k3r << k3v << std::endl;
    //std::cout << k4r << k4v << std::endl;

    p.v = pi->v + 1./6*(k1v + 2.*k2v + 2.*k3v + k4v);
    p.r = pi->r + 1./6*(k1r + 2.*k2r + 2.*k3r + k4r);
    
    //std::cout << pi->r << pi->v << pi->a << std::endl;
    //std::cout << Ei << Bi << std::endl;
    //getchar();
    return p;
}

ParticleData PPSolver::runge4u_mistake(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be, k1r, k2r, k3r, k4r, k1u, k2u, k3u, k4u;
    std::vector<ParticleData>::reverse_iterator pi = (*i)->p.rbegin();
    ParticleData p;
    Vec3<double> pu = momentum_v(pi->v);
    (this->*field_interaction)(Ei, Bi, i);
    
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
    p.r = pi->r + 1./6*(k1r + 2.*k2r + 2.*k3r + k4r);
    return p;
}

ParticleData PPSolver::adams_bashfort4_v(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator p1 = (*i)->p.rbegin(), p2, p3, p4;
    std::vector<ParticleData>::reverse_iterator rend = (*i)->p.rend();
    ParticleData p;
    (this->*field_interaction)(Ei, Bi, i);
    field_extern(Ee, Be, p1->r, t);
    p1->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, p1->v), p1->v);
    //std::cout << p1->r << p1->v << p1->a << std::endl;
    if (p1 != rend) {
        p2 = p1 + 1;
        if (p2 != rend) {
            p3 = p2 + 1;
            if (p3 != rend) {
                p4 = p3 + 1;
                if (p4 != rend) {
                    p.v = p1->v + cdt*(55./24.*p1->a - 59./24.*p2->a + 37./24.*p3->a - 9./24.*p4->a);
                    p.r = p1->r + cdt*(55./24.*p1->v - 59./24.*p2->v + 37./24.*p3->v - 9./24.*p4->v);
                } else {
                    p.v = p1->v + cdt*(23./12.*p1->a - 16./12.*p2->a + 5./12.*p3->a);
                    p.r = p1->r + cdt*(23./12.*p1->v - 16./12.*p2->v + 5./12.*p3->v);
                }
            } else {
                p.v = p1->v + cdt*(1.5*p1->a - 0.5*p2->a);
                p.r = p1->r + cdt*(1.5*p1->v - 0.5*p2->v);
            }
        } else {
            p.v = p1->v + cdt*p1->a;
            p.r = p1->r + cdt*p1->v;
            //std::cout << p.r << p.v << p.a << std::endl;
        }
    }
    //std::cout << p.r << p.v << p.a << std::endl;
    //getchar();
    return p;
}

ParticleData PPSolver::adams_bashfort4_u(std::vector<Particle*>::iterator &i) {
    Vec3<double> Ei, Bi, Ee, Be;
    std::vector<ParticleData>::reverse_iterator p1 = (*i)->p.rbegin(), p2, p3, p4;
    std::vector<ParticleData>::reverse_iterator rend = (*i)->p.rend();
    ParticleData p;
    (this->*field_interaction)(Ei, Bi, i);
    field_extern(Ee, Be, p1->r, t);
    p1->a = acceleration((*i)->type->forcePm(Ee + Ei, Be + Bi, p1->v), p1->v);
    if (p1 != rend) {
        p2 = p1 + 1;
        if (p2 != rend) {
            p3 = p2 + 1;
            if (p3 != rend) {
                p4 = p3 + 1;
                if (p4 != rend) {
                    p.v = velocity_u(momentum_v(p1->v) 
                        + cdt*(55./24.*forcepm(p1->a, p1->v) 
                        - 59./24.*forcepm(p2->a, p2->v) 
                        + 37./24.*forcepm(p3->a, p3->v) 
                        -  9./24.*forcepm(p4->a, p4->v)));
                    p.r = p1->r + cdt*(55./24.*p1->v - 59./24.*p2->v + 37./24.*p3->v - 9./24.*p4->v);
                } else {
                    p.v = velocity_u(momentum_v(p1->v) 
                        + cdt*(23./12.*forcepm(p1->a, p1->v) 
                        - 16./12.*forcepm(p2->a, p2->v) 
                        +  5./12.*forcepm(p3->a, p3->v)));
                    p.r = p1->r + cdt*(23./12.*p1->v - 16./12.*p2->v + 5./12.*p3->v);
                }
            } else {
                p.v = velocity_u(momentum_v(p1->v) 
                    + cdt*(1.5*forcepm(p1->a, p1->v) 
                    - 0.5*forcepm(p2->a, p2->v)));
                p.r = p1->r + cdt*(1.5*p1->v - 0.5*p2->v);
            }
        } else {
            p.v = velocity_u(momentum_v(p1->v) + cdt*forcepm(p1->a, p1->v));
            p.r = p1->r + cdt*p1->v;
        }
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
    size_t result = 0;
    for (std::vector<Injector*>::iterator i = injectors.begin(); i != injectors.end(); i++) {
        //std::cout << "injectors.size() = " << injectors.size() << std::endl;
        (*i)->inject(particles, cdt, nt);
    }
    
    //std::cout << "We are in function" << std::endl;
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
                result++;
                (*i)->nmax = nt;
                (*i)->p.push_back(p);
            }
        }
    } else {

    }
    return result;
}

size_t PPSolver::calcOMP() {
    nt++;
    t += dt;
    ct += cdt;
    size_t result = 0;
    for (std::vector<Injector*>::iterator i = injectors.begin(); i != injectors.end(); i++) {
        //std::cout << "injectors.size() = " << injectors.size() << std::endl;
        (*i)->inject(particles, cdt, nt);
    }
    
    //std::cout << "We are in function" << std::endl;
    if (number_method < RUNGE4V) {
        #pragma omp parallel for schedule(guided) shared(result)
        for (int s = 0; s<(int) particles.size(); s++) {
            std::vector<Particle*>::iterator i = particles.begin() + s;
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
                result++;
                (*i)->nmax = nt;
                (*i)->p.push_back(p);
            }
        }
    } else {

    }
    return result;
}

size_t PPSolver::save(const std::string &filename) {
    //std::cout << "0" << std::endl;
    std::cout << "0" << std::endl;
    if (fd.is_open()) {
        fd.close();
    }
    fd.open(filename.c_str(), std::ofstream::out);
    if (fd.is_open()) {
        //std::cout << "01" << std::endl;
        int s = 0;
        for (std::vector<Particle*>::iterator i = particles.begin(); i != particles.end(); i++) {
            //std::cout << "1" << std::endl;
            fd << "# " << (*i)->nmax << " " << (*i)->p.size() << " " << (*i)->type->type << " ";
            //std::cout << "2" << std::endl;
            (*i)->type->save(fd);
            //std::cout << "3" << std::endl;
            fd << std::endl;
            int k = 0;
            for (std::vector<ParticleData>::iterator j = (*i)->p.begin(); j != (*i)->p.end(); j++) {
                j->save(fd);
                fd << "\n";
                k++;
            }
            s++;
            //std::cout << "Particle number: " << s << " n step: " << k << std::endl;
            fd << std::endl;
        }
        fd.close();
    } else {
        std::cout << "error opening file" << std::endl;
    }
    return 0;
}

size_t PPSolver::load(const std::string &filename) {
    // std::ifstream fd(filename);
    // char c;
    // size_t type;
    // size_t ntype;
    // size_t size;
    // ParticleData p;
    // Particle *ps;
    // ParticleType *ptype;
    // while (fd.eof() == 0) {
    //     fd >> c;
    //     if (c == '#') {
    //         fd >> ntype;
    //         while (ntype > 0) {
    //             fd >> c;

    //             ntype--;
    //         }
    //     }
    //     ps = new Particle();
    //     if (c == '#') {
    //         fd >> ps->nmax >> size >> type;
    //         if (type == 0) {
    //             ptype = new ParticleTypePoint();
    //             ptype->load(fd);
    //         }
    //         if (type == 1) {
    //             ptype = new ParticleTypeBall();
    //             ptype->load(fd);
    //         }
    //         std::cout << "nmax = " << ps->nmax << " size = " << size << std::endl;
    //         for (size_t i = 0; i<size; i++) {
    //             fd >> p.r.x >> p.r.y >> p.r.z >> 
    //                  p.v.x >> p.v.y >> p.v.z >>
    //                  p.a.x >> p.a.y >> p.a.z;
    //             ps->
    //         }
    //     }
    // }
    // fd.close();
    return 0;
}

size_t PPSolver::save(const std::string &filename, const size_t &n) {
    //std::ofstream fd(filename);
    fd.open(filename, std::ofstream::out | std::ofstream::app);
    for (std::vector<Particle*>::iterator i = particles.begin(); i != particles.end(); i++) {
        if (((*i)->nmax >= n) && (n > (*i)->nmax - (*i)->p.size())) {
            size_t index = n - (*i)->nmax + (*i)->p.size() - 1;
            std::vector<ParticleData>::reference p = (*i)->p.at(index);
            fd << n*dt << " " << p.r.x << " " << p.r.y << " " << p.r.z << 
                         " " << p.v.x << " " << p.v.y << " " << p.v.z <<
                         " " << p.a.x << " " << p.a.y << " " << p.a.z << "\n";
        }
    }
    fd << std::endl;
    fd.close();
    return n;
}

void PPSolver::free() {
    for (std::vector<ParticleType*>::iterator i = particle_types.begin(); i != particle_types.end(); i++) {
        delete *i;
    }
    for (std::vector<Injector*>::iterator i = injectors.begin(); i != injectors.end(); i++) {
        delete *i;
    }
    for (std::vector<ExternField*>::iterator i = extern_fields.begin(); i != extern_fields.end(); i++) {
        delete *i;
    }
    for (std::vector<Space3d*>::iterator i = spaces.begin(); i != spaces.end(); i++) {
        delete *i;
    }
    for (std::vector<Particle*>::iterator i = particles.begin(); i != particles.end(); i++) {
        delete *i;
    }
}