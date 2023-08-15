//
// Created by yirao zhang on 2023-08-08.
//

#ifndef COMPOSITEILM_HAMILTONIANMC_H
#define COMPOSITEILM_HAMILTONIANMC_H
#include "ILM.h"
#include "Inference.h"
#include <random>
#include <iostream>
#include <cppad/cppad.hpp>
#endif //COMPOSITEILM_HAMILTONIANMC_H
vector<double> diff_log_target(Individual pop[], double a0, double a1, double beta);
void leapfrog(Individual pop[], vector<double>& z, vector<double>& r, double e, int L);
void HMC(Individual pop[], double samples[][3]);