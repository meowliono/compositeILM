//
// Created by yirao zhang on 2023-06-25.
//

#ifndef COMPOSITEILM_MCMC_H
#define COMPOSITEILM_MCMC_H

#endif //COMPOSITEILM_MCMC_H
#include "ILM.h"
#include "Inference.h"
#include <random>
#include <iostream>
#include <boost/math/distributions/gamma.hpp>
void MH(Individual pop[], int t_end, double samples[][3]);
void MH_log(Individual pop[], int t_end, double samples[][3]);