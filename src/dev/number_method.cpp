#include "number_method.h"

#ifdef HAS_TEST

class _Data {
public:
    double x;
    double y;

    _Data() : x(0), y(0) {};
    _Data(const double &x, const double &y) : x(x), y(y) {};

    _Data operator +(const _Data &a) const {
        return _Data(x + a.x, y + a.y);
    }

    _Data operator *(const double &a) const {
        return _Data(x*a, y*a);
    }

    static _Data dxdt(const _Data &a, const double &t) {
        return _Data(-a.y, a.x);
    }

    static inline size_t size() { return 2; }
};

void test_number_method() {
    double dt = 3e-3;
    size_t N = 10000;
    
    std::vector<NumberMethod<_Data>*> methods;
    NumberMethod<_Data>* method;
    method = new MethodEuler<_Data>(dt, _Data::dxdt);
    methods.push_back(method);
    method = new MethodRunge4<_Data>(dt, _Data::dxdt);
    methods.push_back(method);
    std::cout << "Number of methods = " << methods.size() << std::endl;
    
    std::ofstream file = std::ofstream("results.gp", std::ofstream::out);
    file << "plot cos(x) w l lw 2 title 'Theory'";
    for (size_t i = 0; i<methods.size(); i++) {
        std::cout << methods[i]->getName() << std::endl;
        file << ",\\\n     'results.dat' " << "using 1:" << i*_Data::size() + 2 << " w l lw 2 title '" << methods[i]->getName() << "'";
    }
    file.close();
    file = std::ofstream("results.dat", std::ofstream::out);
    file << std::scientific;
    _Data* u = new _Data[methods.size()];
    file << "\n" << 0.;
    for (size_t i = 0; i<methods.size(); i++) {
        u[i].x = 1;
        u[i].y = 0;
        file << " " << u[i].x << " " << u[i].y;
    } 
    for (size_t i = 0; i<N; i++) {
        file << "\n" << (i + 1)*dt;
        for (size_t j = 0; j<methods.size(); j++) {
            u[j] = methods[j]->next(u[j], i*dt);
            file << " " << u[j].x << " " << u[j].y;
        } 
    }

    file.close();
}

#endif // HAS_TEST