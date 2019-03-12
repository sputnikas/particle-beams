#pragma once

#include <iostream>
#include <string>

// Data has 2 operator Data + Data, and double * Data
template <typename Data>
class NumberMethod {
    std::string name;
public:
    NumberMethod(std::string name) : name(name) {};

    std::string getName() { return name; }

    virtual Data next(const Data &data, const double &t) = 0;
};

template <typename Data>
class MethodEuler : public NumberMethod<Data> {  
    Data (*dxdt)(const Data &data, const double &t);
    double dt;
public:
    MethodEuler() : NumberMethod<Data>("Euler"), dt(1) { };
    MethodEuler(double dt) : NumberMethod<Data>("Euler"), dt(dt) { };
    MethodEuler(double dt, Data (*dxdt)(const Data &, const double &)) : NumberMethod<Data>("Euler"), dxdt(dxdt), dt(dt)  { };
    
    Data next(const Data &data, const double &t) {
        return data + dxdt(data, t)*dt;
    }
};

template <typename Data>
class MethodRunge4 : public NumberMethod<Data> {
    Data (*dxdt)(const Data &data, const double &t);
    double dt;
public:
    MethodRunge4() : NumberMethod<Data>("Runge4"), dt(1) { };
    MethodRunge4(double dt) : NumberMethod<Data>("Runge4"), dt(dt)  { };
    MethodRunge4(double dt, Data (*dxdt)(const Data &, const double &)) : NumberMethod<Data>("Runge4"), dxdt(dxdt), dt(dt) { };
    
    Data next(const Data &data, const double &t) {
        Data k1, k2, k3, k4;
        k1 = dxdt(data, t)*dt;
        k2 = dxdt(data + k1*0.5, t + dt/2)*dt;
        k3 = dxdt(data + k2*0.5, t + dt/2)*dt;
        k4 = dxdt(data + k3, t + dt)*dt;
        return data + (k1 + k2*2. + k3*2. + k4)*(1./6);
    }
};

#ifdef HAS_TEST

#include <fstream>
#include <vector>

void test_number_method();

#endif // HAS_TEST