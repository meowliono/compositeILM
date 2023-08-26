//
// Created by yirao zhang on 2023-08-09.
//

#ifndef COMPOSITEILM_MCMC_COMPOSITE_H
#define COMPOSITEILM_MCMC_COMPOSITE_H
#include "ILM.h"
#include "Inference_composite.h"
#include <random>
#include <iostream>
#include <boost/math/distributions/gamma.hpp>
#endif //COMPOSITEILM_MCMC_COMPOSITE_H
void MH_composite(map<int, vector<Individual>>clusters, double samples[][4]);
void MH_composite_K(map<int, vector<Individual>>clusters, double samples[][6]);
void MH_composite_tv(map<int, vector<Individual>>clusters, double samples[][5]);
void MH_composite_bc(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, double samples[][5]);