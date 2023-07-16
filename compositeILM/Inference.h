//
// Created by yirao zhang on 2023-07-12.
//

#ifndef COMPOSITEILM_INFERENCE_H
#define COMPOSITEILM_INFERENCE_H

#endif //COMPOSITEILM_INFERENCE_H
#include "ILM.h"
extern double prior_param[];
extern double log_lhd;
void set_prior_Gaussian(double param[6]);
double calculate_log_likelihood(Individual pop[], double a0, double a1, double n1);
double post_dist_dens(Individual pop[],double a0, double a1, double n1);