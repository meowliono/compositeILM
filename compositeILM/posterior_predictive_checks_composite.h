//
// Created by yirao zhang on 2023-10-14.
//

#ifndef COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_COMPOSITE_H
#define COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_COMPOSITE_H
#include "posterior_predictive_checks.h"
#include "Clustering.h"
#include <vector>
#endif //COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_COMPOSITE_H
double prob_sus_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, int id, double s_a0, double s_a1, double s_beta, double t_beta);
double prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta);
double prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta);
void status_change_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, int id, double s_a0, double s_a1, double s_beta, double t_beta);
void status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta);
void status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta);
void pop_status_change_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, double s_a0, double s_a1, double s_beta, double t_beta);
void pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, double s_a0, double s_a1, double s_beta, double t_beta);
void pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta);
void simulate_ppc_afterstartday_bc_tilde_woeps(Individual* pop, vector<DataPoint> centroids, int start_day, int end_day);
void simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual* pop, int start_day, int end_day);
void simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual* pop, int start_day, int end_day);
void partially_simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual* pop, int step_size);;
void partially_simulate_ppc_afterstartday_bc_tilde_woeps(Individual* pop, vector<DataPoint> centroids, const vector<int>& records);