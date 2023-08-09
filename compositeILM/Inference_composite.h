//
// Created by yirao zhang on 2023-08-09.
//

#ifndef COMPOSITEILM_INFERENCE_COMPOSITE_H
#define COMPOSITEILM_INFERENCE_COMPOSITE_H
#include "ILM.h"
#include "Inference.h"
#include <boost/math/distributions/gamma.hpp>
#include <map>
#include <vector>
#endif //COMPOSITEILM_INFERENCE_COMPOSITE_H
using namespace std;
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps);
double composite_log_likelihood(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps);
double composite_post_dens(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps);