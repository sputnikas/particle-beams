#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "vec3.h"

const double EQ   = 1.6021892E-19;
const double EM   = 9.109534E-31;
const double EPS0 = 8.854187817E-12;
const double K0   = 8.987551788E9;
const double CL   = 2.99792458E8;
const double RL   = 2.8179402894E-15;
const double CL2  = 8.9875517873681764E16;

//////////////////////////////////////////////////////////////////////////////
// Функции для СТО
// скорости в единицах c, energy - в единицах энергии покоя, momentum - в единицах mc
//////////////////////////////////////////////////////////////////////////////

inline double energy_v(const Vec3<double> &v) {
    return 1/sqrt(1 - v.norm2());
}

inline Vec3<double> momentum_v(const Vec3<double> &v) {
    return v/sqrt(1 - v.norm2());
}

inline double energy_u(const Vec3<double> &u) {
    return sqrt(1 + u.norm2());
}

inline Vec3<double> velocity_u(const Vec3<double> &u) {
    return u/sqrt(1 + u.norm2());
}

inline Vec3<double> acceleration(const Vec3<double> &fpm, const Vec3<double> &v) {
    return (fpm - v*fpm.dot(v))*sqrt(1 - v.norm2());
}

inline Vec3<double> forcepm(const Vec3<double> &acceleration, const Vec3<double> &v) {
    return (acceleration + acceleration.dot(v)*v)/sqrt(1 - v.norm2());
}

//////////////////////////////////////////////////////////////////////////////
// Particle
// После долгих размышлений появилась одна идея:
//      частицы, должны быть упорядочены во временной массив!
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// ParticleData
// хранит информацию о положении частицы скоростях и ускорении
// здесь снова находится краегульный камень
// дело в том, что если рассматривать частицы, как нечто обобщённое,
// то они бывают различными:
//      по скорости и характеру взаимодействия:
//          - нерелятивистские
//          - релятивисткие
//      по типу распределения заряда, массы, формы и соответственно поля их окружающего:
//          - с простой симметрией (r, v, a) и всё
//              - точечные
//              - шаровые
//              - сферические
//          - с более сложной симметрией (появляются дополнительные степени свободы, как углы Эйлера или сложнее)
//              - точечные диполи
//              - точечные квадруполи ...
//              - кольца
//              - диски
//              - стержни
//          - с дополнительными необычными степенями свободы:
//              - частицы переменной массы (они ж крупные, могут и разделиться при столкновении!)
// и как учитывать всё это многообразие непонятно...
// а ведь ещё и можно задачу решать как в переменных (r, v, a), так и в переменных (r, p, F)!
// и в данном месте я окончательно выпадаю в осадок:
//      чтобы сделать такую абстракцию, которую потом легко расширить, изменить - нужно быть гением c++
//      единственный вариант, который я пока вижу для данной задачи - это использование специализированных шаблонов
//      и да мне такая сложность не нравится
// поэтому мы пока остановимся на одном единственном типе частиц с точечной (шаровой) симметрией
//////////////////////////////////////////////////////////////////////////////

class ParticleData {
public:
    Vec3<double> r; // в [м]
    Vec3<double> v; // в единицах скорости света
    Vec3<double> a; // в [м^-1], то есть отнесённое к скорости света в квадрате

    size_t load(std::istream &f);
    size_t save(std::ostream &f);
    size_t load(char *f);
    size_t save(FILE *f);
};

//////////////////////////////////////////////////////////////////////////////
// ParticleType
// нужен чтобы хранить параметры общие для групп частиц
// поля в field прибавляются к текущему значению
// E и B измеряются в [В/м], то есть вместо B рассматривается произведение cB
// A и phi измеряются в [В], то есть вместо A рассматривается Ac
//////////////////////////////////////////////////////////////////////////////

enum PARTICLE_TYPES {
    PARTICLE_TYPE_POINT,
    PARTICLE_TYPE_BALL,
    PARTICLE_TYPE_SPHERE,
    PARTICLE_TYPE_NUMBER
};

class ParticleType {
public:
    size_t type;
    ParticleType(const size_t &type);

    virtual ~ParticleType();
    virtual void field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) = 0;
    virtual void potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) = 0;
    virtual void fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r) = 0;
    virtual void potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r) = 0;
    virtual Vec3<double> forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v) = 0;
    virtual size_t load(std::istream &f) = 0;
    virtual size_t save(std::ostream &f) = 0;
    virtual size_t load(char *f) = 0;
    virtual size_t save(FILE *f) = 0;
};

class ParticleTypePoint : public ParticleType {
public:
    double q;
    double m;
    double qpm;

    ParticleTypePoint();
    ParticleTypePoint(const double &q, const double &m);

    void field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    void fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    Vec3<double> forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v);
    size_t load(std::istream &f);
    size_t save(std::ostream &f);
    size_t load(char *f);
    size_t save(FILE *f);
};

class ParticleTypeSphere : public ParticleType {
public:
    double q;
    double m;
    double qpm;
    double radius;

    ParticleTypeSphere();
    ParticleTypeSphere(const double &q, const double &m, const double &radius);

    void field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    void fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    Vec3<double> forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v);
    size_t load(std::istream &f);
    size_t save(std::ostream &f);
    size_t load(char *f);
    size_t save(FILE *f);
};

class ParticleTypeBall : public ParticleType {
public:
    double q;
    double m;
    double qpm;
    double radius;

    ParticleTypeBall();
    ParticleTypeBall(const double &q, const double &m, const double &radius);

    void field(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potential(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    void fieldN(Vec3<double> &E, Vec3<double> &B, const ParticleData &p, const Vec3<double> &r);
    void potentialN(Vec3<double> &A, double &phi, const ParticleData &p, const Vec3<double> &r);
    Vec3<double> forcePm(const Vec3<double> &E, const Vec3<double> &B, const Vec3<double> &v);
    size_t load(std::istream &f);
    size_t save(std::ostream &f);
    size_t load(char *f);
    size_t save(FILE *f);
};

//////////////////////////////////////////////////////////////////////////////
// Particle - in time
// частицы расположены в прямом порядке:
// last в момент времени nmax
// 1 в момент времени nmax - size() - то есть в начальный момент времени
//////////////////////////////////////////////////////////////////////////////

class Particle {
public:
    std::vector<ParticleData> p;
    size_t nmax;
    ParticleType *type;
};

///////////////////////////////////////////////////////////////////////////////
// Sphere method
// Различные реализации метода сфер
// - метод интервалов
// - метод бисекции
// - метод бисекции с уточнением
///////////////////////////////////////////////////////////////////////////////

std::pair<size_t, ParticleData> sm_left_intervals(Particle *p, const Vec3<double> &r, const double &t, const double &dt);
std::pair<size_t, ParticleData> sm_left_bisection(Particle *p, const Vec3<double> &r, const double &t, const double &dt);
std::pair<size_t, ParticleData> sm_bisection(Particle *p, const Vec3<double> &r, const double &t, const double &dt);
// Тестирование метода сфер
void test_sphere_method();
