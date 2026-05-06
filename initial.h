// initial.h
#pragma once
#include <vector>

double initial_wavepacket(double x, int m);
double exact_solution(double x, double t, int m);
// 新增：谱初值及精确解
double initial_spectrum(double x, const std::vector<double>& psi);
double exact_solution_spectrum(double x, double t, const std::vector<double>& psi);