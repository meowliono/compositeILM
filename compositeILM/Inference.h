//
// Created by yirao zhang on 2023-07-12.
//

#ifndef COMPOSITEILM_INFERENCE_H
#define COMPOSITEILM_INFERENCE_H

#endif //COMPOSITEILM_INFERENCE_H
#include "ILM.h"
#include <boost/math/distributions/gamma.hpp>
#include <map>
#include <vector>
using namespace std;
extern vector<double> prior_param;
void set_prior_parameters(double param[], int num);
double calculate_log_likelihood(Individual pop[], int t_end, double a0, double a1, double beta);
double calculate_log_likelihood_logprior(Individual pop[], int t_end, double ta0, double ta1, double tbeta);
double post_dist_dens(Individual pop[],int t_end, double a0, double a1, double beta);
double post_dist_dens_logprior(Individual pop[],int t_end, double ta0, double ta1, double tbeta);