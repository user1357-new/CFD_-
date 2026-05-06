// initial.cpp
#include "initial.h"
#include "constants.h"
#include <cmath>
#include <vector>

double initial_wavepacket(double x, int m) {
    double sum = 0.0;
    for (int l = 1; l <= m; ++l)
        sum += std::sin(2.0 * PI * l * x);
    return sum / m;
}

double exact_solution(double x, double t, int m) {
    return initial_wavepacket(x - A * t, m);
}
double initial_spectrum(double x, const std::vector<double>& psi) {
    const double eps = 0.1;
    const int Kmax = 64;
    const double k0 = 24.0;
    double sum = 0.0;
    for (int k = 1; k <= Kmax; ++k) {
        double Ek = std::pow(k / k0, 4.0) * std::exp(-2.0 * (k / k0) * (k / k0));
        sum += std::sqrt(Ek) * std::sin(2.0 * PI * k * (x + psi[k - 1]));
    }
    return 1.0 + eps * sum;
}

double exact_solution_spectrum(double x, double t, const std::vector<double>& psi) {
    return initial_spectrum(x - A * t, psi);
}