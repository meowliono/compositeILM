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
extern double prior_param[];
extern double log_lhd;
void set_prior_parameters(double param[]);
double calculate_log_likelihood(Individual pop[], double a0, double a1, double beta);
double calculate_log_likelihood_logprior(Individual pop[], double ta0, double ta1, double tbeta);
double post_dist_dens(Individual pop[],double a0, double a1, double beta);
double post_dist_dens_logprior(Individual pop[],double ta0, double ta1, double tbeta);
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps);
double composite_log_likelihood(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps);
double composite_post_dens(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps);