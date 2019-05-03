#pragma once

#include <cmath>
#include <cstdlib>

inline double arsh(double x) {
    return log(x*x + sqrt(x*x + 1));
}

inline double arch(double x) {
    return log(x*x + sqrt(x*x - 1));
}

inline double arth(double x) {
    return 0.5*log((x + 1)/(x - 1));
}

// Возвращает равномерно распр. случайное число от 0 до 1
inline double random() {
    return ((double) rand())/RAND_MAX;
}