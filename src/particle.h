#include <iostream>
#include <vector>

#include "vec3.h"
/*
class ParticleSliceBase { 
};

template <class ParticleData>
class ParticleSlice : ParticleSliceBase {
    std::vector<ParticleData> data;
};

class ParticleDataBase {
    int type;
public:
    ParticleDataBase() : type(0) {};
    ParticleDataBase(const int &type) : type(type) {};

    int getType() const { return type; };
    void setType(const int &t) { type = t; };
};

class ParticlePointR : public ParticleDataBase {
    static std::vector<double> m;
    static std::vector<double> q;
    static std::vector<double> qdm;

    Vec3<double> r;
    Vec3<double> v;
    Vec3<double> a;

    ParticlePointR() 
        : ParticleDataBase(), r(Vec3<double>()), v(Vec3<double>()), a(Vec3<double>()) {};
    ParticlePointR(const int &type) 
        : ParticleDataBase(type), r(Vec3<double>()), v(Vec3<double>()), a(Vec3<double>()) {};
    ParticlePointR(const Vec3<double> &r) 
        : ParticleDataBase(), r(r), v(Vec3<double>()), a(Vec3<double>()) {};
    ParticlePointR(const int &type, const Vec3<double> &r) 
        : ParticleDataBase(type), r(r), v(Vec3<double>()), a(Vec3<double>()) {};
    ParticlePointR(const Vec3<double> &r, const Vec3<double> &v) 
        : ParticleDataBase(), r(r), v(v), a(Vec3<double>()) {};
    ParticlePointR(const int &type, const Vec3<double> &r, const Vec3<double> &v) 
        : ParticleDataBase(type), r(r), v(v), a(Vec3<double>()) {};
    ParticlePointR(const ParticlePointR &p) 
        : ParticleDataBase(p.getType()), r(p.r), v(p.v), a(p.a) {};
    ~ParticlePointR() {};

    void interaction(Vec3<double> &E, Vec3<double> &B, Vec3<double> rs) {
        Vec3<double> R = rs - r;
        double Rn = norm(R);
        Vec3<double> n = R/Rn;
        Vec3<double> nvc = n - v/CL;
        E = q[getType()]*K0*(nvc*(1 - norm2(v)/CL2) + Rn/CL2*n*(nvc*a))/Rn/Rn/pow(1 - dot(n,v)/CL);
        B = E*n/CL;
    }   
};

class ParticlePoint : ParticleDataBase {
    Vec3<double> r;
    Vec3<double> v;
};

class ParticlePoint : ParticleDataBase {
    Vec3<double> r;
    Vec3<double> p;
};

*/

const double EQ = 1.6021892E-19;
const double EM = 9.109534E-31;
const double EPS0 = 8.854187817E-12;
const double K0 = 8.987551788E9;
const double CL = 2.99792458E8;
const double RL = 2.8179402894E-15;
const double CL2 = 8.9875517873681764E16;

//////////////////////////////////////////////////////////////////////////////
// Particle
// пока 1 тип - точечные или шаровые частицы с точно определённым типом взаимодействия
// по Лиенару-Вихерту, закону Кулона и Био-Савара-Лапласа, с редукцией на малых
// расстояниях 
//////////////////////////////////////////////////////////////////////////////

class Particle {
public:
    int id;
    Vec3<double> r; // в [м]
    Vec3<double> v; // в c
    Vec3<double> a; // 
};

typedef std::vector<Particle>::iterator ParticleIndex;

class ParticleSlice {
public:
    static std::vector<int> nmin;
    static std::vector<int> nmax;
    static std::vector<std::vector<int>> id;
    static std::vector<int> type;
    std::vector<Particle> particles;

    inline int size() { return particles.size(); }
    inline void add(const Particle &p) {
        particles.push_back(p);
    }
    inline void rm(const ParticleIndex &i) {
        particles.erase(i);
    }
    inline ParticleIndex begin() {
        return particles.begin();
    }
    inline ParticleIndex end() {
        return particles.end();
    }
    static inline int all_particles() { return nmin.size(); }
};

class ParticleType {
public:
    virtual void field(Vec3<double> &E, Vec3<double> &B, const Particle &p, const Vec3<double> &r) = 0;
};

class ParticlePointR : ParticleType {
public:
    double m;
    double q;
    double qdm;
    
    void field(Vec3<double> &E, Vec3<double> &B, const Particle &p, const Vec3<double> &r) {
        Vec3<double> R = r - p.r;
        Vec3<double> beta = p.v / CL;
        double Rn = norm(R);
        Vec3<double> n = R/Rn;
        Vec3<double> nvc = n - beta;
        E = q*K0*(nvc*(1 - norm2(beta)) + Rn/CL2*n*(nvc*p.a))/Rn/Rn/pow(1 - dot(n,beta), 1.5);
        B = E*n/CL;
    }
};

class ParticlePointN : ParticleType {
public:
    double m;
    double q;
    double qdm;
    
    void field(Vec3<double> &E, Vec3<double> &B, const Particle &p, const Vec3<double> &r) {
        Vec3<double> R = r - p.r;
        Vec3<double> beta = p.v / CL;
        double Rn = norm(R);
        Vec3<double> n = R/Rn;
        Vec3<double> nvc = n - beta;
        E = q*K0*(nvc*(1 - norm2(beta)) + Rn/CL2*n*(nvc*p.a))/Rn/Rn/pow(1 - dot(n,beta), 1.5);
        B = E*n/CL;
    }
};

//////////////////////////////////////////////////////////////////////////////
//  Injector
//////////////////////////////////////////////////////////////////////////////

#include <sstream>

class Injector {
public:
    double I;   // ток протекающий через инжектор
    double rho; // плотность пространственного заряда
    double kgec;   // коэффициент укрупнения
    int n;      // количество инжектируемых частиц за шаг
    int type;   // тип частицы
    double tmin;
    double tmax;

    virtual void inject(ParticleSlice *particles, const int &step, const double &dt) = 0;
    virtual double area() = 0;
};

inline double random(const double &xmin, const double &xmax){
    return xmin + (xmax - xmin)*1.*rand()/RAND_MAX;
}

class InjectorRectangularZUniform : public Injector{
public:
    double wx;
    double wy;
    double x;
    double y;
    double z;
    double vz;

    void inject(ParticleSlice *slice, const int &step, const double &dt) {
        double t = step*dt;
        if (((t < tmax) && (t > tmin)) || (tmin*tmax == 0)) {
            for (int i = 0; i<n; i++) {
                Particle p;
                p.r.x = random(x - wx/2, x + wx/2);
                p.r.y = random(y - wy/2, y + wy/2);
                p.r.z = random(z, z + vz*dt);
                p.v.x = 0;
                p.v.y = 0;
                p.v.z = vz;
                p.a.x = 0;
                p.a.y = 0;
                p.a.z = 0;
                p.id = ParticleSlice::all_particles();
                slice->add(p);
                ParticleSlice::nmin.push_back(step);
                std::vector<int> id;
                id.push_back(slice->size());
                ParticleSlice::id.push_back(id);
                slice->type.push_back(type);
            }
        }
    };
    double area() { return wx*wy; };
};

std::ostream &operator<< (std::ostream &stream, const InjectorRectangularZUniform &inj) {
    stream << "inj.wx = " << inj.wx << "\n";
    stream << "inj.wy = " << inj.wy << "\n";
    stream << "inj.x = "  << inj.x  << "\n";
    stream << "inj.y = "  << inj.y  << "\n";
    stream << "inj.z = "  << inj.z  << "\n";
    stream << "inj.vz = " << inj.vz << "\n";
    stream << "inj.I = "  << inj.I  << "\n";
    stream << "inj.rho = "   << inj.rho   << "\n"; 
    stream << "inj.kgec = "  << inj.kgec  << "\n";
    stream << "inj.n = "  << inj.n  << "\n";
    stream << "inj.tmin = " << inj.tmin << "\n";
    stream << "inj.tmax = " << inj.tmax << "\n";
}

//////////////////////////////////////////////////////////////////////////////
//  Space
//////////////////////////////////////////////////////////////////////////////

class Space {
public:
    virtual int in(const Vec3<double> &r) = 0;
};

class SpaceCube {
public:
    Vec3<double> min;
    Vec3<double> max;

    int in(const Vec3<double> &r) {
        return ((r.x > min.x) && (r.y > min.y) && (r.z > min.z) && (r.x < max.x) && (r.y < max.y) && (r.z < max.z));
    }
}; 

class SpaceAbovePlane {
public:
    Vec3<double> point;
    Vec3<double> normal;

    int in(const Vec3<double> &r) {
        return (dot(r - point, normal) > 0);
    }
};

class SpaceBetweenPlanes {
public:
    Vec3<double> point;
    Vec3<double> normal;
    double h;

    int in(const Vec3<double> &r) {
        double d = dot(r - point, normal);
        return ((d > 0) && (d < h));
    }
};

//////////////////////////////////////////////////////////////////////////////
//  ExternField
//////////////////////////////////////////////////////////////////////////////

class ExternField {
public:
    Space *space;
    virtual void field(Vec3<double> &E, Vec3<double> &B, const Vec3<double> &r, const double &t) = 0;
};

class ExternFieldConst : public ExternField {
public:
    Vec3<double> E0;
    Vec3<double> B0;

    void field(Vec3<double> &E, Vec3<double> &B, const Vec3<double> &r, const double &t) {
        if (space->in(r)) {
            E = E0;
            B = B0;
        } else {
            E = Vec3<double>();
            B = Vec3<double>();
        }
    };
};

namespace Relativity {
    inline double energy(const Vec3<double> &u) {
        return sqrt(norm2(u) + 1);
    }

    inline Vec3<double> vel(const Vec3<double> &u) {
        return u/sqrt(norm2(u) + 1);
    }

    inline Vec3<double> momentum(const Vec3<double> &v) {
        return v/sqrt(1 - norm2(v));
    }
}

//////////////////////////////////////////////////////////////////////////////
//  PPSolver
//////////////////////////////////////////////////////////////////////////////

class PPSolver {
public:
    typedef std::vector<Injector*>::iterator      InjectorIndex;
    typedef std::vector<Space*>::iterator         SpaceIndex;
    typedef std::vector<ExternField*>::iterator   ExternFieldIndex;

    std::vector<ParticleType*>  types;
    std::vector<Injector*>      injectors;
    std::vector<Space*>         spaces;
    std::vector<ExternField*>   efield;
    std::vector<ParticleSlice*> slices;
    
    double dt;
    int step;

    int (*sm)(Particle &p, const Vec3<double> &r, const int &pnumber, const double &t, const double &dt);
    int (*nm)(Particle *p, const int &step, const double &dt, const Vec3<double> &Ei, const Vec3<double> &Bi, const double &t);

    void euler(Particle *p, const int &step, const double &dt, const Vec3<double> &Ei, const Vec3<double> &Bi, const double &t) {
        p->v = Relativity::vel(Relativity::momentum(p->v) + eforces(Ei, Bi, p->r, p->v, t)*dt);
        p->r += p->v*dt; 
    }

    int runge4(Particle *p, const int &step, const double &dt, const Vec3<double> &Ei, const Vec3<double> &Bi, const double &t) {
        Vec3<double> rk1, rk2, rk3, rk4, pk1, pk2, pk3, pk4;
        rk1 = Relativity::vel(Relativity::momentum(p->v) + eforces(Ei, Bi, p->r, p->v, t)*dt);
    }

    inline Particle *getParticle(const int &nt, const int &pnumber) {
        return &slices[nt]->particles[ParticleSlice::id[pnumber][nt - ParticleSlice::nmin[pnumber]]];
    }

    int sm_intervals(Particle &p, const Vec3<double> &r, const int &pnumber, const double &t, const double &dt) {
        int nmax = ParticleSlice::nmax[pnumber];
        int step = (int) floor(t/dt);
        if (nmax == -1) nmax = step;
        int nmin = ParticleSlice::nmin[pnumber];
        if (norm<double>(r - getParticle(nmin, pnumber)->r) >= (t - nmin*dt)*CL)
            return -1;
        if (norm<double>(r - getParticle(nmin, pnumber)->r) < (t - (nmax + 1)*dt)*CL)
            return -1;
        int n_sphere = step - nmax;
        int n_max = nmax;
        while (n_max >= nmin){
            if (norm<double>(r - getParticle(nmax, pnumber)->r) < n_sphere*CL*dt)
                return n_max;
            else {
                n_max--;
                n_sphere++;
            }
        }
        return -1;
    }

    void interaction(Vec3<double> &E, Vec3<double> &B, const Vec3<double> &r) {
        E = Vec3<double>();
        B = Vec3<double>();
        for (int j = 0; j<slices[step]->all_particles(); j++) {
            Vec3<double> Ei = Vec3<double>();
            Vec3<double> Bi = Vec3<double>();
            Particle p;
            if (sm(p, r, j, step, dt) != -1) {
                types[ParticleSlice::type[p.id]]->field(Ei, Bi, p, r);
            };
            E += Ei;
            B += Bi;
        }
    }

    Vec3<double> eforces(const Vec3<double> &Ei, const Vec3<double> &Bi, const Vec3<double> &r, const Vec3<double> &v, const double &t) {
        Vec3<double> E = Ei;
        Vec3<double> B = Bi;
        for (ExternFieldIndex j = efield.begin(); j != efield.end(); j++) {
            Vec3<double> Ee = Vec3<double>();
            Vec3<double> Be = Vec3<double>();
            (*j)->field(Ee, Be, r, step*dt);
            E += Ee;
            B += Be;
        }
        return E + v*B;
    }

    void next() {
        double t = step*dt;
        ParticleSlice *slice = new ParticleSlice();
        for (ParticleIndex i = slices[step]->begin(); i != slices[step]->end(); i++) {
            slice->add(*i);
        }
        for (InjectorIndex i = injectors.begin(); i != injectors.end(); i++) {
            (*i)->inject(slice, step*dt, dt);
        }
        slices.push_back(slice);
        Vec3<double> E, B;
        for (ParticleIndex i = slices[step]->begin(); i != slices[step]->end(); i++) {
            
        }
        step++;
    }

}