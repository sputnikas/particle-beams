#include "particle2.h"

size_t ParticleData::load(std::istream &f) {
    f >> r.x >> r.y >> r.z >> v.x >> v.y >> v.z >> a.x >> a.y >> a.z;
    return 0;
}

size_t ParticleData::save(std::ostream &f) {
    f << r.x << " " << r.y << " " << r.z << " " <<
         v.x << " " << v.y << " " << v.z << " " <<
         a.x << " " << a.y << " " << a.z;
    return 0;
}

size_t ParticleData::load(char *f) {
    sscanf(f, "%le %le %le %le %le %le %le %le %le", &r.x, &r.y, &r.z, &v.x, &v.y, &v.z, &a.x, &a.y, &a.z);
    return 0;
}

size_t ParticleData::save(FILE *f) {
    fprintf(f, "%le %le %le %le %le %le %le %le %le", r.x, r.y, r.z, v.x, v.y, v.z, a.x, a.y, a.z);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
// ParticleType
// нужен чтобы хранить параметры общие свойства для групп частиц
// поля в field прибавляются к текущему значению
//////////////////////////////////////////////////////////////////////////////

ParticleType::ParticleType(const size_t &type) : type(type) {}
ParticleType::~ParticleType() {}

//////////////////////////////////////////////////////////////////////////////
// relativistic point type
//////////////////////////////////////////////////////////////////////////////

ParticleTypePoint::ParticleTypePoint() : ParticleType(0), q(EQ), m(EM), qpm(EQ/EM) {}
ParticleTypePoint::ParticleTypePoint(const double &q, const double &m) : ParticleType(PARTICLE_TYPE_POINT), q(q), m(m), qpm(q/m) {}

void ParticleTypePoint::field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn != 0.0) {
        Vec3<double> n = R/Rn;
        double p1 = 1 - n.dot(p.v);
        double k = K0*q/Rn/Rn/p1/p1/p1;
        Vec3<double> Ep = k*((n - p.v)*(1 - p.v.norm2() + R.dot(p.a)) - Rn*p.a*p1);
        E += Ep;
        B += n.cross(Ep);
    }
}

void ParticleTypePoint::potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn != 0.0) {
        double p1 = Rn - R.dot(p.v);
        double Ap = K0*q/p1;
        A += Ap*p.v;
        phi += Ap;
    }
}

void ParticleTypePoint::fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn != 0.0) {
        Vec3<double> n = R/Rn;
        Vec3<double> Ep = K0*q/Rn/Rn*n;
        E += Ep;
        B += p.v.cross(Ep);
    }
}

void ParticleTypePoint::potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn != 0.0) {
        double phip = K0*q/Rn;
        A += phip*p.v;
        phi +=  phip;
    }
}

Vec3<double> ParticleTypePoint::forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v) {
    return qpm*(E + v.cross(B))/CL2;
}

size_t ParticleTypePoint::load(std::istream &f) {
    f >> q >> m >> qpm;
    return 0;
}

size_t ParticleTypePoint::save(std::ostream &f) {
    f << q << " " << m << " " << qpm;
    return 0;
}

size_t ParticleTypePoint::load(char *f) {
    sscanf(f, "%le %le %le", &q, &m, &qpm);
    return 0;
}

size_t ParticleTypePoint::save(FILE *f) {
    fprintf(f, "%le %le %le", q, m, qpm);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
// relativistic ball type
//////////////////////////////////////////////////////////////////////////////

ParticleTypeBall::ParticleTypeBall() : ParticleType(1), q(EQ), m(EM), qpm(EQ/EM), radius(RL) {}
ParticleTypeBall::ParticleTypeBall(const double &q, const double &m, const double &radius) : ParticleType(PARTICLE_TYPE_BALL), q(q), m(m), qpm(q/m), radius(radius) {}

void ParticleTypeBall::field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    double Rss = Rn - R.dot(p.v);
    if (Rss >= radius*sqrt(1 - p.v.norm2())) {
        Vec3<double> n = R/Rn;
        double k = K0*q/Rss/Rss/Rss*Rn;
        Vec3<double> Ep = k*((n - p.v)*(1 - p.v.norm2() + R.dot(p.a)) - p.a*Rss);
        E += Ep;
        B += n.cross(Ep);
    } else {
        //std::cout << "Rss = " << Rss << std::endl;
        double k = K0*q/radius/Rss;
        double k1 = 1.0 - p.v.norm2() + R.dot(p.a);
        double gamma2 = 1/(1 - p.v.norm2());
        double gamma = sqrt(gamma2);
        double gamma3 = gamma*gamma2;
        double x = Rss/radius;
        double phi = (1.5 - 0.5*x*x*gamma2)*gamma;
        double dphi_dxi = -Rss/radius/radius*gamma3;
        double dphi_dct = -1.5*gamma3*(1.0 - x*x*gamma2)*p.v.dot(p.a);
        Vec3<double> Ep = k*((- k1*dphi_dxi + dphi_dct)*(R - Rn*p.v) - phi*p.a*Rn);
        E += Ep;
        B += R.cross(Ep)/Rn;
    }
}

void ParticleTypeBall::potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    double Rss = Rn - R.dot(p.v);
    if (Rss >= radius*sqrt(1 - p.v.norm2())) {
        double Ap = K0*q/Rss;
        A += Ap*p.v;
        phi += Ap;
    } else {
        double g = 1.0 - p.v.norm2();
        double Ap = K0*q*(1.5 - 0.5*Rss*Rss/radius/radius/g)/radius/sqrt(g);
        A += Ap*p.v;
        phi += Ap;
    }
}

void ParticleTypeBall::fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn >= radius) {
        Vec3<double> Ep = K0*q/Rn/Rn/Rn*R;
        E += Ep;
        B += p.v.cross(Ep);
    } else {
        Vec3<double> Ep = K0*q/radius/radius/radius*R;
        E += Ep;
        B += p.v.cross(Ep);
    }
}

void ParticleTypeBall::potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) {
    Vec3<double> R = r - p.r;
    double Rn = R.norm();
    if (Rn >= radius) {
        double phip = K0*q/Rn;
        A += phip*p.v;
        phi +=  phip;
    } else {
        double Ap = K0*q*(1.5 - 0.5*Rn*Rn/radius/radius)/radius;
        A += Ap*p.v;
        phi += Ap;
    }
}

Vec3<double> ParticleTypeBall::forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v) {
    return qpm*(E + v.cross(B))/CL2;
}

size_t ParticleTypeBall::load(std::istream &f) {
    f >> q >> m >> qpm >> radius;
    return 0;
}

size_t ParticleTypeBall::save(std::ostream &f) {
    f << q << " " << m << " " << qpm << " " << radius;
    return 0;
}

size_t ParticleTypeBall::load(char *f) {
    sscanf(f, "%le %le %le %le", &q, &m, &qpm, &radius);
    return 0;
}

size_t ParticleTypeBall::save(FILE *f) {
    fprintf(f, "%le %le %le %le", q, m, qpm, radius);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Sphere method
// Различные реализации метода сфер
// - метод интервалов
// - метод бисекции
// - метод бисекции с уточнением
///////////////////////////////////////////////////////////////////////////////

std::pair<size_t, ParticleData> sm_left_intervals(  Particle *p,
                                                    const Vec3<double> &r,
                                                    const double &t,
                                                    const double &dt    )
{
    ParticleData result;
    std::vector<ParticleData>::iterator pl = p->p.begin();
    std::vector<ParticleData>::reverse_iterator pr = p->p.rbegin();
    double ct(CL*t);
    double cdt(CL*dt);
    size_t nmin(p->nmax - p->p.size() + 1);
    size_t nmax(p->nmax);
    double f;

    //std::cout << "nmin = " << nmin << " nmax = " << nmax << std::endl;
    //std::cout << "end : "    << p->p.end()->r    << std::endl;
    //std::cout << "begin : "  << p->p.begin()->r  << std::endl;
    //std::cout << "rend : "   << p->p.rend()->r   << std::endl;
    //std::cout << "rbegin : " << p->p.rbegin()->r << std::endl;

    if ((f = (ct - nmin*cdt) - (r - pl->r).norm()) < 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    //std::cout << "n = " << nmin << " f = " << f << std::endl;
    if ((f = (ct - nmax*cdt) - (r - pr->r).norm()) > 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    //std::cout << "n = " << nmax << " f = " << f << std::endl;
    do {
        f = (ct - nmax*cdt) - (r - pr->r).norm();
        nmax--;
        pr++;
        //getchar();
    } while (pr != (p->p.rend()) && (f < 0));
    //std::cout << "n = " << nmax + 1 << " f = " << f << std::endl;
    return std::pair<size_t, ParticleData>(nmax + 1, p->p[nmax + 1 - nmin]);
}

std::pair<size_t, ParticleData> sm_left_bisection(  Particle *p,
                                                    const Vec3<double> &r,
                                                    const double &t,
                                                    const double &dt    )
{
    ParticleData result;
    std::vector<ParticleData>::iterator pl = p->p.begin();
    std::vector<ParticleData>::reverse_iterator pr = p->p.rbegin();
    double ct(CL*t);
    double cdt(CL*dt);
    size_t nmin(p->nmax - p->p.size() + 1);
    size_t nmax(p->nmax);
    size_t n;
    double f;

    if ((f = (ct - nmax*cdt) - (r - pr->r).norm()) > 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    if ((f = (ct - nmin*cdt) - (r - pl->r).norm()) < 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    size_t n2min = nmin;
    do {
        n = ((nmin + nmax) % 2 == 0) ? (nmin + nmax)/2 : nmin/2 + nmax/2;
        f = (ct - n*cdt) - (r - p->p[n - n2min].r).norm();
        if (f == 0) {
            return std::pair<size_t, ParticleData>(n + n2min, p->p[n - n2min]);
        }
        if (f > 0) {
            nmin = n;
        } else {
            nmax = n;
        }
    } while (nmax != nmin + 1);
    return std::pair<size_t, ParticleData>(n2min + nmin, p->p[nmin - n2min]);
}

std::pair<size_t, ParticleData> sm_bisection(   Particle *p,
                                                const Vec3<double> &r,
                                                const double &t,
                                                const double &dt    )
{
    ParticleData result;
    std::vector<ParticleData>::iterator pl = p->p.begin();
    std::vector<ParticleData>::reverse_iterator pr = p->p.rbegin();
    double ct(CL*t);
    double cdt(CL*dt);
    size_t nmin(p->nmax - p->p.size() + 1);
    size_t nmax(p->nmax);
    size_t n;
    double f;

    if ((f = (ct - nmax*cdt) - (r - pr->r).norm()) > 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    if ((f = (ct - nmin*cdt) - (r - pl->r).norm()) < 0) {
        return std::pair<int, ParticleData>(-1, result);
    }
    size_t n2min = nmin;
    do {
        n = ((nmin + nmax) % 2 == 0) ? (nmin + nmax)/2 : nmin/2 + nmax/2;
        f = (ct - n*cdt) - (r - p->p[n - n2min].r).norm();
        if (f == 0) {
            return std::pair<size_t, ParticleData>(n + n2min, p->p[n - n2min]);
        }
        if (f > 0) {
            nmin = n;
        } else {
            nmax = n;
        }
        //std::cout << "nmin = " << nmin << "nmax = " << nmax << " f = " << f << std::endl;
        //getchar();
    } while (nmax != nmin + 1);

    ParticleData *p1(&(p->p[nmin - n2min])), *p2(&(p->p[nmin + 1 - n2min]));
    double cdt2 =  cdt*cdt;
    double cdt3 = cdt2*cdt;
    double cdt4 = cdt3*cdt;
    double cdt5 = cdt4*cdt;
    const double e = 0.01;
    Vec3<double> A3(  10.*(p2->r - p1->r)/cdt3 - (4.*p2->v + 6.*p1->v)/cdt2 + (0.5*p2->a - 1.5*p1->a)/cdt);
    Vec3<double> A4(- 15.*(p2->r - p1->r)/cdt4 + (7.*p2->v + 8.*p1->v)/cdt3 - (    p2->a - 1.5*p1->a)/cdt2);
    Vec3<double> A5(   6.*(p2->r - p1->r)/cdt5 - 3.*(p2->v +    p1->v)/cdt4 + 0.5*(p2->a -     p1->a)/cdt3);
    double dl = 0, dr = 1, d, d2, d3, d4, d5;
    do {
        d = (dl + dr)/2;
        d2 = d*d;
        d3 = d2*d;
        d4 = d3*d;
        d5 = d4*d;
        result.r = p1->r + p1->v*d*cdt +  p1->a*d2*cdt2 +     A3*d3*cdt3 +    A4*d4*cdt4 + A5*d5*cdt5;
        f = (ct - (nmin + d)*cdt) - (r - result.r).norm();
        if (f == 0) {
            break;
        }
        if (f > 0) {
            dl = d;
        } else {
            dr = d;
        }
    } while (dr - dl > e);
    result.v = p1->v + p1->a*d*cdt +  3.*A3*d2*cdt2 +  4.*A4*d3*cdt3 + 5.*A5*d4*cdt4;
    result.a = p1->a + 6.*A3*d*cdt + 12.*A4*d2*cdt2 + 20.*A5*d3*cdt3;

    if (result.v.norm2() > 1) {
        std::cout << "\nnmin = " << nmin << ": " << p1->r << p1->v << p1->a << std::endl;
        std::cout << "nmax = " << nmax << ": " << p2->r << p2->v << p2->a << std::endl;
        std::cout << result.r << result.v << result.a << std::endl;
        getchar();
    }
    return std::pair<size_t, ParticleData>(n2min + nmin, result);
}

///////////////////////////////////////////////////////////////////////////////
// Тестирование метода сфер
///////////////////////////////////////////////////////////////////////////////

void test_sphere_method_l(Particle *p, const Vec3<double> &r, const double &t, const double &dt, const Vec3<double> &v0, const Vec3<double> &r0) {
    // Теоретические значения промежутка времени
    Vec3<double> a = (r - r0)/CL;
    double b = (t - a.dot(v0))/(1 - v0.norm2());
    double c = (t*t - a.norm2())/(1 - v0.norm2());
    double t1 = b - sqrt(b*b - c);
    double t2 = b + sqrt(b*b - c);
    std::cout << "t1 = " << t1 << "; t = " << t << "; t2 = " << t2 << std::endl;
    std::cout << "n1 = " << t1/dt << "; n = " << t/dt << "; n2 = " << t2/dt  << std::endl;
    std::cout << "0: t = " << t << "; n = " << t1/dt <<  "; r = " << r0 + v0*t1*CL << "; v = " << v0 << "; a = " << Vec3<double>() << std::endl;

    // Результат работы методов сфер
    std::pair<size_t, ParticleData> pr;
    pr = sm_left_intervals(p, r, t, dt);
    std::cout << "1: t = " << t << "; n = " << pr.first <<  "; r = " << pr.second.r << "; v = " << pr.second.v << "; a = " << pr.second.a << std::endl;
    //getchar();
    pr = sm_left_bisection(p, r, t, dt);
    std::cout << "2: t = " << t << "; n = " << pr.first <<  "; r = " << pr.second.r << "; v = " << pr.second.v << "; a = " << pr.second.a << std::endl;
    //getchar();
        if (pr.first != (size_t) -1) {
        Vec3<double> E, B;
        E = Vec3<double>();
        B = Vec3<double>();
        p->type->field(E, B, pr.second, r);
        std::cout << "delay rel: E = " << E << " B = " << B << std::endl;

        E = Vec3<double>();
        B = Vec3<double>();
        p->type->fieldN(E, B, pr.second, r);
        std::cout << "delay non: E = " << E << " B = " << B << std::endl;

        E = Vec3<double>();
        B = Vec3<double>();
        p->type->fieldN(E, B, p->p[0], r);
        std::cout << "non delay: E = " << E << " B = " << B << std::endl;

        Vec3<double> A;
        double phi;
        A = Vec3<double>();
        phi = 0;
        p->type->potential(A, phi, pr.second, r);
        std::cout << "delay rel: A = " << A << " phi = " << phi << std::endl;

        A = Vec3<double>();
        phi = 0;
        p->type->potentialN(A, phi, pr.second, r);
        std::cout << "delay non: A = " << A << " phi = " << phi << std::endl;

        A = Vec3<double>();
        phi = 0;
        p->type->potentialN(A, phi, p->p[0], r);
        std::cout << "non delay: A = " << A << " phi = " << phi << std::endl;
    } else {
        std::cout << "Particle for field out of time range!" << std::endl;
    }
    pr = sm_bisection(p, r, t, dt);
    std::cout << "3: t = " << t << "; n = " << pr.first <<  "; r = " << pr.second.r << "; v = " << pr.second.v << "; a = " << pr.second.a << std::endl;
    //getchar();

    if (pr.first != (size_t) -1) {
        Vec3<double> E, B;
        E = Vec3<double>();
        B = Vec3<double>();
        p->type->field(E, B, pr.second, r);
        std::cout << "delay rel: E = " << E << " B = " << B << std::endl;

        E = Vec3<double>();
        B = Vec3<double>();
        p->type->fieldN(E, B, pr.second, r);
        std::cout << "delay non: E = " << E << " B = " << B << std::endl;

        E = Vec3<double>();
        B = Vec3<double>();
        p->type->fieldN(E, B, p->p[0], r);
        std::cout << "non delay: E = " << E << " B = " << B << std::endl;

        Vec3<double> A;
        double phi;
        A = Vec3<double>();
        phi = 0;
        p->type->potential(A, phi, pr.second, r);
        std::cout << "delay rel: A = " << A << " phi = " << phi << std::endl;

        A = Vec3<double>();
        phi = 0;
        p->type->potentialN(A, phi, pr.second, r);
        std::cout << "delay non: A = " << A << " phi = " << phi << std::endl;

        A = Vec3<double>();
        phi = 0;
        p->type->potentialN(A, phi, p->p[0], r);
        std::cout << "non delay: A = " << A << " phi = " << phi << std::endl;
    } else {
        std::cout << "Particle for field out of time range!" << std::endl;
    }
    getchar();
}

void test_sphere_method() {
    int n = 1000;
    double q = 1000000*EQ;
    double m = 1000000*EM;
    Vec3<double> v0 = Vec3<double>(0.0, 0.0, 0.1);
    Vec3<double> r0 = Vec3<double>(0.0, 0.0, 0.0);
    double dt = 1e-12;

    double cdt = CL*dt;
    Particle p;
    ParticleData pd;
    double ct = 0;
    p.nmax = 0;
    for (int i = 0; i<n; i++) {
        pd.r = r0 + v0*ct;
        pd.v = v0;
        ct += cdt;
        p.p.push_back(pd);
        p.nmax = i;
        //std::cout << "r = " << pd.r << std::endl;
    }
    std::cout << "nmax = " << p.nmax << std::endl;

    Vec3<double> r1 = Vec3<double>( 0.0, 0.01, 0.0);
    Vec3<double> r2 = Vec3<double>(0.01,  0.0, 0.0);
    Vec3<double> r3 = Vec3<double>(0.01, 0.01, 0.0);
    double t1 = 10*dt;
    double t2 = 100*dt;
    double t3 = 200*dt;
    double t4 = 2000*dt;

    p.type = new ParticleTypePoint(q, m);
    test_sphere_method_l(&p, r1, t1, dt, v0, r0);
    test_sphere_method_l(&p, r1, t2, dt, v0, r0);
    test_sphere_method_l(&p, r1, t3, dt, v0, r0);
    test_sphere_method_l(&p, r1, t4, dt, v0, r0);
    delete p.type;

    p.type = new ParticleTypeBall(q, m, 0.1);
    test_sphere_method_l(&p, r1, t1, dt, v0, r0);
    test_sphere_method_l(&p, r1, t2, dt, v0, r0);
    test_sphere_method_l(&p, r1, t3, dt, v0, r0);
    test_sphere_method_l(&p, r1, t4, dt, v0, r0);
    delete p.type;

    p.p.clear();
}

