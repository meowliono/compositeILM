//
// Created by yirao zhang on 2023-08-09.
//

#ifndef COMPOSITEILM_INFERENCE_COMPOSITE_H
#define COMPOSITEILM_INFERENCE_COMPOSITE_H
#include "ILM.h"
#include "Inference.h"
#include "Clustering.h"
#include <boost/math/distributions/gamma.hpp>
#include <map>
#include <vector>
#endif //COMPOSITEILM_INFERENCE_COMPOSITE_H
using namespace std;
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps);
double cluster_log_likelihood_tv(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps, double delta);
double cluster_log_likelihood_bc(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int id, double a0, double a1, double beta, double eps, double tilde_beta);
double composite_log_likelihood(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps);
double composite_log_likelihood(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps[],bool ClusterDepedent);
double composite_log_likelihood_tv(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps, double delta);
double composite_log_likelihood_bc(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, double a0, double a1, double beta, double eps, double tilde_beta);
double composite_post_dens(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps);
double composite_post_dens(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps[], bool ClusterDependent);
double composite_post_dens_tv(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps, double delta);
double composite_post_dens_bc(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, double a0, double a1, double beta, double eps, double tilde_beta);