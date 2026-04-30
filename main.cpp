// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "constants.h"
#include "initial.h"
#include "schemes.h"

// ---------- 线性格式的通用模拟函数 ----------
double run_linear_sim(const Grid& grid, int m,
                      void (*rhs_func)(std::vector<double>&, const std::vector<double>&, const Grid&)) {
    const double dt_max = CFL * grid.dx;
    std::vector<double> u(grid.N);
    for (int i = 0; i < grid.N; ++i)
        u[i] = initial_wavepacket(grid.x[i], m);

    double t = 0.0;
    RK4 rk4;
    while (t < T_END - 1e-12) {
        double dt = std::min(dt_max, T_END - t);
        rk4.step(u, grid, dt, rhs_func);
        t += dt;
    }
    return compute_L2_error(u, grid, T_END, m);
}

// ---------- SA-DRP 专用模拟函数 ----------
double run_sadrp_sim(const Grid& grid, int m) {
    const double dt_max = CFL * grid.dx;
    std::vector<double> u(grid.N);
    for (int i = 0; i < grid.N; ++i)
        u[i] = initial_wavepacket(grid.x[i], m);

    double t = 0.0;
    while (t < T_END - 1e-12) {
        double dt = std::min(dt_max, T_END - t);
        step_SADRP(u, grid, dt);
        t += dt;
    }
    return compute_L2_error(u, grid, T_END, m);
}

int main() {
    std::ofstream ferr("error_convergence.csv");
    ferr << "Scheme,m,N,L2_error\n";

    std::vector<int> Nvals = { 64, 128, 256, 512, 1024 };
    std::vector<int> mvals = {20};

    // 定义要测试的格式（只保留 DRP 和 MDCD，DRP-M 系数不确切可省略）
    std::vector<std::pair<std::string, void(*)(std::vector<double>&, const std::vector<double>&, const Grid&)>> schemes = {
        {"DRP",  rhs_DRP},
        {"DRP-M", rhs_DRPM},
        {"MDCD", rhs_MDCD}
    };

    for (auto& s : schemes) {
        for (int m : mvals) {
            for (int N : Nvals) {
                Grid grid(N);
                double err = run_linear_sim(grid, m, s.second);
                ferr << s.first << "," << m << "," << N << "," << err << "\n";
                std::cout << s.first << " m=" << m << " N=" << N << " L2=" << err << std::endl;
            }
        }
    }

    // SA-DRP 单独运行
    for (int m : mvals) {
        for (int N : Nvals) {
            Grid grid(N);
            double err = run_sadrp_sim(grid, m);
            ferr << "SA-DRP," << m << "," << N << "," << err << "\n";
            std::cout << "SA-DRP m=" << m << " N=" << N << " L2=" << err << std::endl;
        }
    }
    ferr.close();

    // ---------- 输出波形对比 (N=256, m=20) ----------
    const int Nw = 256;
    const int mw = 20;
    Grid gw(Nw);
    std::ofstream fwave("waveform_N256_m20.csv");
    fwave << "x,Exact,DRP,MDCD,SA-DRP\n";

    // 计算各格式的最终波形
    auto compute_wave = [&](auto rhs_func) {
        const double dt_max = CFL * gw.dx;
        std::vector<double> u(Nw);
        for (int i = 0; i < Nw; ++i) u[i] = initial_wavepacket(gw.x[i], mw);
        double t = 0.0;
        RK4 rk4;
        while (t < T_END - 1e-12) {
            double dt = std::min(dt_max, T_END - t);
            rk4.step(u, gw, dt, rhs_func);
            t += dt;
        }
        return u;
    };

    auto u_drp  = compute_wave(rhs_DRP);
    auto u_mdcd = compute_wave(rhs_MDCD);

    // SA-DRP 波形
    std::vector<double> u_sadrp(Nw);
    for (int i = 0; i < Nw; ++i) u_sadrp[i] = initial_wavepacket(gw.x[i], mw);
    double t = 0.0;
    double dt_max = CFL * gw.dx;
    while (t < T_END - 1e-12) {
        double dt = std::min(dt_max, T_END - t);
        step_SADRP(u_sadrp, gw, dt);
        t += dt;
    }

    for (int i = 0; i < Nw; ++i) {
        fwave << gw.x[i] << ","
              << exact_solution(gw.x[i], T_END, mw) << ","
              << u_drp[i] << "," << u_mdcd[i] << "," << u_sadrp[i] << "\n";
    }
    fwave.close();

    std::cout << "\nFiles generated: error_convergence.csv, waveform_N256_m20.csv\n";
    return 0;
}