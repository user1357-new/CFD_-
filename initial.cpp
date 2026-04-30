// initial.cpp
#include "initial.h"
#include "constants.h"
#include <cmath>

double initial_wavepacket(double x, int m) {
    double sum = 0.0;
    for (int l = 1; l <= m; ++l)
        sum += std::sin(2.0 * PI * l * x);
    return sum / m;
}

double exact_solution(double x, double t, int m) {
    return initial_wavepacket(x - A * t, m);
}