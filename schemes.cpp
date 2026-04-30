// schemes.cpp
#include "schemes.h"
#include "initial.h"
#include <cmath>
#include <algorithm>

// ================= Grid =================
Grid::Grid(int n) : N(n), dx(1.0 / n), x(n) {
    for (int i = 0; i < N; ++i) x[i] = i * dx;
}
int Grid::idx(int i) const {
    return (i % N + N) % N;
}

// ================= RK4 =================
void RK4::step(std::vector<double>& u, const Grid& grid,
               double dt,
               void (*rhs)(std::vector<double>&, const std::vector<double>&, const Grid&)) {
    std::vector<double> k1(u.size()), k2(u.size()), k3(u.size()), k4(u.size());
    std::vector<double> utmp(u.size());

    rhs(k1, u, grid);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + 0.5 * dt * k1[i];

    rhs(k2, utmp, grid);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + 0.5 * dt * k2[i];

    rhs(k3, utmp, grid);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + dt * k3[i];

    rhs(k4, utmp, grid);
    for (size_t i = 0; i < u.size(); ++i)
        u[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}
// ================= DRP =================
const double DRP_C[4] = {0.0, 0.79926643, -0.18941314, 0.02651995};

void rhs_DRP(std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid) {
    const int N = grid.N;
    const double idx = 1.0 / grid.dx;
    for (int i = 0; i < N; ++i) {
        double d = 0.0;
        for (int j = 1; j <= 3; ++j)
            d += DRP_C[j] * (u[grid.idx(i + j)] - u[grid.idx(i - j)]);
        dudt[i] = -A * d * idx;
    }
}

// ================= DRP-M =================
const double DRPM_C[4] = {0.0, 0.806394, -0.196077, 0.032244};

void rhs_DRPM(std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid) {
    const int N = grid.N;
    const double idx = 1.0 / grid.dx;
    for (int i = 0; i < N; ++i) {
        double d = 0.0;
        for (int j = 1; j <= 3; ++j)
            d += DRPM_C[j] * (u[grid.idx(i + j)] - u[grid.idx(i - j)]);
        dudt[i] = -A * d * idx;
    }
}

// ================= MDCD =================
const double GAMMA_DISP_MDCD = 0.0463783;
const double GAMMA_DISS_MDCD = 0.001;

void rhs_MDCD(std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid) {
    const int N = grid.N;
    const double dx = grid.dx;
    const double idx = 1.0 / dx;
    for (int i = 0; i < N; ++i) {
        auto hat = [&](int j) {
            int ja  = grid.idx(j);
            int jb  = grid.idx(j+1);
            int jm1 = grid.idx(j-1);
            int jp2 = grid.idx(j+2);
            int jm2 = grid.idx(j-2);
            int jp3 = grid.idx(j+3);
            double r = 0.0;
            r += (0.5 * GAMMA_DISP_MDCD + 0.5 * GAMMA_DISS_MDCD) * u[jm2];
            r += (-1.5 * GAMMA_DISP_MDCD - 2.5 * GAMMA_DISS_MDCD - 1.0/12.0) * u[jm1];
            r += (GAMMA_DISP_MDCD + 5.0 * GAMMA_DISS_MDCD + 7.0/12.0) * u[ja];
            r += (GAMMA_DISP_MDCD - 5.0 * GAMMA_DISS_MDCD + 7.0/12.0) * u[jb];
            r += (-1.5 * GAMMA_DISP_MDCD + 2.5 * GAMMA_DISS_MDCD - 1.0/12.0) * u[jp2];
            r += (0.5 * GAMMA_DISP_MDCD - 0.5 * GAMMA_DISS_MDCD) * u[jp3];
            return r * idx;
        };
        double fp = hat(i);      // j+1/2
        double fm = hat(i-1);    // j-1/2
        dudt[i] = -(fp - fm);
    }
}

// ================= SA-DRP =================
void compute_interface_params(const std::vector<double>& u, const Grid& grid,
                              std::vector<double>& gdisp, std::vector<double>& gdiss) {
    const int N = grid.N;
    gdisp.resize(N);
    gdiss.resize(N);
    for (int i = 0; i < N; ++i) { // i 对应界面 j+1/2 (j = i)
        int im2 = grid.idx(i-2), im1 = grid.idx(i-1), i0 = i;
        int ip1 = grid.idx(i+1), ip2 = grid.idx(i+2), ip3 = grid.idx(i+3);

        double umax = u[im2], umin = u[im2];
        for (auto idx : {im1, i0, ip1, ip2, ip3}) {
            if (u[idx] > umax) umax = u[idx];
            if (u[idx] < umin) umin = u[idx];
        }
        double kESW = 0.0;
        double gamma_disp = 1.0/30.0;
        double gamma_diss = 0.001;

        if (umax - umin > 1e-3) {
            double S1 = u[ip1] - 2.0*u[i0] + u[im1];
            double S2 = (u[ip2] - 2.0*u[i0] + u[im2]) / 4.0;
            double S3 = u[ip2] - 2.0*u[ip1] + u[i0];
            double S4 = (u[ip3] - 2.0*u[ip1] + u[im1]) / 4.0;
            double C1 = u[ip1] - u[i0];
            double C2 = (u[ip2] - u[im1]) / 3.0;

            auto term = [&](double a, double b, double eps) {
                double abs_sum = std::fabs(a + b);
                double abs_diff = std::fabs(a - b);
                return (std::fabs(abs_sum - abs_diff) + eps) / (abs_sum + abs_diff + eps);
            };
            double t1 = term(S1, S2, EPS);
            double t2 = term(S3, S4, EPS);
            double t3 = (std::fabs(std::fabs(C1+C2) - 0.5*std::fabs(C1-C2)) + 2.0*EPS)
                       / (std::fabs(C1+C2) + std::fabs(C1-C2) + EPS);
            double val = 2.0 * std::min({t1, t2, t3}) - 1.0;
            if (val > 1.0) val = 1.0;
            if (val < -1.0) val = -1.0;
            kESW = std::acos(val);

            if (kESW < 0.01) {
                gamma_disp = 1.0/30.0;
            } else if (kESW < K_C_DISP) {
                double sk = std::sin(kESW);
                double s2k = std::sin(2.0*kESW);
                double s3k = std::sin(3.0*kESW);
                double num = kESW + (1.0/6.0)*s2k - (4.0/3.0)*sk;
                double den = s3k - 4.0*s2k + 5.0*sk;
                gamma_disp = num / den;
            } else {
                gamma_disp = GAMMA_DISP_CLIP;
            }

            if (kESW <= 1.0) {
                gamma_diss = 0.001;
            } else {
                double tmp = 0.001 + 0.011 * std::sqrt((kESW - 1.0) / (PI - 1.0));
                gamma_diss = std::min(tmp, 0.012);
            }
        }
        gdisp[i] = gamma_disp;
        gdiss[i] = gamma_diss;
    }
}

void rhs_SADRP_cached(std::vector<double>& dudt, const std::vector<double>& u,
                      const Grid& grid,
                      const std::vector<double>& gamma_disp,
                      const std::vector<double>& gamma_diss) {
    const int N = grid.N;
    const double idx = 1.0 / grid.dx;
    for (int i = 0; i < N; ++i) {
        // 通量函数 (局部 lambda)
        auto hat = [&](double gd, double gs, int im2, int im1, int i0, int ip1, int ip2, int ip3) {
            double r = 0.0;
            r += (0.5*gd + 0.5*gs) * u[im2];
            r += (-1.5*gd - 2.5*gs - 1.0/12.0) * u[im1];
            r += (gd + 5.0*gs + 7.0/12.0) * u[i0];
            r += (gd - 5.0*gs + 7.0/12.0) * u[ip1];
            r += (-1.5*gd + 2.5*gs - 1.0/12.0) * u[ip2];
            r += (0.5*gd - 0.5*gs) * u[ip3];
            return r * idx;
        };
        // j+1/2 界面参数 (索引 i)
        double gd_p = gamma_disp[i];
        double gs_p = gamma_diss[i];
        int im2 = grid.idx(i-2), im1 = grid.idx(i-1), i0 = i;
        int ip1 = grid.idx(i+1), ip2 = grid.idx(i+2), ip3 = grid.idx(i+3);
        double fp = hat(gd_p, gs_p, im2, im1, i0, ip1, ip2, ip3);

        // j-1/2 界面参数 (索引 i-1)
        int i_prev = grid.idx(i-1);
        double gd_m = gamma_disp[i_prev];
        double gs_m = gamma_diss[i_prev];
        int jm2 = grid.idx(i-3), jm1 = grid.idx(i-2), j0 = grid.idx(i-1);
        int jp1 = grid.idx(i), jp2 = grid.idx(i+1), jp3 = grid.idx(i+2);
        double fm = hat(gd_m, gs_m, jm2, jm1, j0, jp1, jp2, jp3);

        dudt[i] = -(fp - fm);
    }
}

void step_SADRP(std::vector<double>& u, const Grid& grid, double dt) {
    std::vector<double> gdisp, gdiss;
    compute_interface_params(u, grid, gdisp, gdiss);  // 仅第一阶段计算一次

    std::vector<double> k1(u.size()), k2(u.size()), k3(u.size()), k4(u.size()), utmp(u.size());

    rhs_SADRP_cached(k1, u, grid, gdisp, gdiss);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + 0.5*dt*k1[i];

    rhs_SADRP_cached(k2, utmp, grid, gdisp, gdiss);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + 0.5*dt*k2[i];

    rhs_SADRP_cached(k3, utmp, grid, gdisp, gdiss);
    for (size_t i = 0; i < u.size(); ++i) utmp[i] = u[i] + dt*k3[i];

    rhs_SADRP_cached(k4, utmp, grid, gdisp, gdiss);
    for (size_t i = 0; i < u.size(); ++i)
        u[i] += dt/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}

// ================= 通用误差计算 =================
double compute_L2_error(const std::vector<double>& u, const Grid& grid, double t, int m) {
    double err = 0.0;
    for (int i = 0; i < grid.N; ++i) {
        double diff = u[i] - exact_solution(grid.x[i], t, m);
        err += diff * diff;
    }
    return std::sqrt(err * grid.dx);
}