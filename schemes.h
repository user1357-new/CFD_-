// schemes.h
#pragma once
#include <vector>
#include "constants.h"

class Grid {
public:
    int N;
    double dx;
    std::vector<double> x;
    Grid(int n);
    int idx(int i) const;
};

class RK4 {
public:
    // 修改：步长由外部传入，不再作为成员变量
    void step(std::vector<double>& u, const Grid& grid,
              double dt,
              void (*rhs)(std::vector<double>&, const std::vector<double>&, const Grid&));
};

void rhs_DRP (std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid);
void rhs_DRPM(std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid);
void rhs_MDCD(std::vector<double>& dudt, const std::vector<double>& u, const Grid& grid);

void compute_interface_params(const std::vector<double>& u, const Grid& grid,
                              std::vector<double>& gamma_disp,
                              std::vector<double>& gamma_diss);
void rhs_SADRP_cached(std::vector<double>& dudt, const std::vector<double>& u,
                      const Grid& grid,
                      const std::vector<double>& gamma_disp,
                      const std::vector<double>& gamma_diss);
void step_SADRP(std::vector<double>& u, const Grid& grid, double dt);

double compute_L2_error(const std::vector<double>& u, const Grid& grid, double t, int m);