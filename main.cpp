// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include "constants.h"
#include "initial.h"
#include "schemes.h"

int main() {
    // 固定随机种子，保证可复现
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<double> psi(64);
    for (int k = 0; k < 64; ++k)
        psi[k] = dist(rng);

    // 准备收敛数据文件
    std::ofstream fconv("error_convergence.csv");
    fconv << "N,DRP,DRP-M,MDCD,SA-DRP,UPWIND1,UPWIND2,UPWIND3\n";
    //fconv << "N,UPWIND1,UPWIND2,UPWIND3\n";
    // 存放N=256时的波形解（用于最后写文件）
    std::vector<double> u_drp_256, u_drpm_256, u_mdcd_256, u_sadrp_256;
    //std::vector<double> u_upwind1_256, u_upwind2_256, u_upwind3_256;

    // 循环不同网格点数
    std::vector<int> Nvals = {2048,4096};
    // 用于计算收敛阶：存储上一个N的误差（七种格式）
    double prev_err[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //double prev_err[4] = {0.0, 0.0, 0.0, 0.0};
    int prev_N = 0;

    for (int N : Nvals) {
        Grid grid(N);
        const double dt_max = CFL * grid.dx;

        // 通用线性格式模拟函数 (DRP, DRP-M, MDCD, UPWIND1, UPWIND2, UPWIND3)
        auto run_linear = [&](void (*rhs)(std::vector<double>&, const std::vector<double>&, const Grid&),
                              std::vector<double>& u_out) {
            std::vector<double> u(N);
            for (int i = 0; i < N; ++i)
                u[i] = initial_spectrum(grid.x[i], psi);

            double t = 0.0;
            RK4 rk4;
            while (t < T_END - 1e-12) {
                double dt = std::min(dt_max, T_END - t);
                rk4.step(u, grid, dt, rhs);
                t += dt;
            }
            u_out = u;
            return compute_L1_error(u, grid, T_END, psi);
        };

        // DRP
        std::vector<double> u_drp;
        double err_drp = run_linear(rhs_DRP, u_drp);

        // DRP-M
        std::vector<double> u_drpm;
        double err_drpm = run_linear(rhs_DRPM, u_drpm);

        // MDCD
        std::vector<double> u_mdcd;
        double err_mdcd = run_linear(rhs_MDCD, u_mdcd);

        // 一阶迎风
        std::vector<double> u_upwind1;
        double err_upwind1 = run_linear(rhs_UPWIND1, u_upwind1);

        // 二阶迎风
        std::vector<double> u_upwind2;
        double err_upwind2 = run_linear(rhs_UPWIND2, u_upwind2);

        // 三阶迎风
        std::vector<double> u_upwind3;
        double err_upwind3 = run_linear(rhs_UPWIND3, u_upwind3);

        // SA-DRP 专用
        std::vector<double> u_sadrp(N);
        for (int i = 0; i < N; ++i)
            u_sadrp[i] = initial_spectrum(grid.x[i], psi);
        double t = 0.0;
        while (t < T_END - 1e-12) {
            double dt = std::min(dt_max, T_END - t);
            step_SADRP(u_sadrp, grid, dt);
            t += dt;
        }
        double err_sadrp = compute_L1_error(u_sadrp, grid, T_END, psi);

        // 保存至CSV
        fconv << N << "," << err_drp << "," << err_drpm << "," << err_mdcd << "," 
              << err_sadrp << "," << err_upwind1 << "," << err_upwind2 << "," << err_upwind3 << "\n";
       
       

        // 控制台输出当前N误差及收敛阶（N>64时计算）
        std::cout << "N=" << N
                 << " | DRP: " << err_drp
                  << " | DRP-M: " << err_drpm
                  << " | MDCD: " << err_mdcd
                  << " | SA-DRP: " << err_sadrp
                  << " | UPWIND1: " << err_upwind1
                  << " | UPWIND2: " << err_upwind2
                  << " | UPWIND3: " << err_upwind3;
        if (prev_N > 0) {
            double ratio = std::log(2.0);  // N倍半，阶 = log(err1/err2)/log(2)
            std::cout << " | order (vs prev N): "
                      << std::log(prev_err[0] / err_drp) / ratio << ", "
                      << std::log(prev_err[1] / err_drpm) / ratio << ", "
                      << std::log(prev_err[2] / err_mdcd) / ratio << ", "
                      << std::log(prev_err[3] / err_sadrp) / ratio << ", "
                      << std::log(prev_err[4] / err_upwind1) / ratio << ", "
                      << std::log(prev_err[5] / err_upwind2) / ratio << ", "
                      << std::log(prev_err[6] / err_upwind3) / ratio;
        }
        std::cout << std::endl;

        // 保存N=256波形数据
        if (N == 256) {
            u_drp_256   = u_drp;
            u_drpm_256  = u_drpm;
            u_mdcd_256  = u_mdcd;
            u_sadrp_256 = u_sadrp;
            u_upwind1_256 = u_upwind1;
            u_upwind2_256 = u_upwind2;
            u_upwind3_256 = u_upwind3;
        }

        // 更新上一个误差记录
        prev_err[0] = err_drp;
        prev_err[1] = err_drpm;
        prev_err[2] = err_mdcd;
        prev_err[3] = err_sadrp;
        prev_err[4] = err_upwind1;
        prev_err[5] = err_upwind2;
        prev_err[6] = err_upwind3;
        prev_N = N;
    }
    fconv.close();

    // 输出 N=256 波形文件
    std::ofstream fwave("waveform_spectrum_N256.csv");
    fwave << "x,Exact,DRP,DRP-M,MDCD,SA-DRP,UPWIND1,UPWIND2,UPWIND3\n";
    // 重建grid(256)仅用于输出坐标（也可直接使用dx=1/256）
    Grid g256(256);
    for (int i = 0; i < 256; ++i) {
        fwave << g256.x[i] << ","
              << exact_solution_spectrum(g256.x[i], T_END, psi) << ","
              << u_drp_256[i] << ","
              << u_drpm_256[i] << ","
              << u_mdcd_256[i] << ","
              << u_sadrp_256[i] << ","
              << u_upwind1_256[i] << ","
              << u_upwind2_256[i] << ","
              << u_upwind3_256[i] << "\n";
    }
    fwave.close();

    std::cout << "\nConvergence data saved to error_convergence.csv\n";
    std::cout << "Waveform (N=256) saved to waveform_spectrum_N256.csv\n";
    return 0;
}