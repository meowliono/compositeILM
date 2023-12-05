//
// Created by yirao zhang on 2023-09-18.
//

#ifndef COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_H
#define COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_H
#include "ILM.h"
#include <random>
#include <map>
#include <vector>
#endif //COMPOSITEILM_POSTERIOR_PREDICTIVE_CHECKS_H
extern int ppc_crt_time;
extern int ppc_inf_time_tab[][300];
extern int ppc_rem_time_tab[][300];
double prob_sus_ppc(Individual population[], int id, double a0, double a1, double beta, double eps);
double prob_sus_ppc_afterstartday(Individual population[], int i, int id, double a0, double a1, double beta, double eps);
void status_change_ppc_beforestartday(Individual population[], int id, double a0, double a1, double beta, double eps);
void status_change_ppc_afterstartday(Individual population[], int i, int id, double a0, double a1, double beta, double eps);
void pop_status_change_ppc_beforestartday(Individual population[], double a0, double a1, double beta, double eps);
void pop_status_change_ppc_afterstartday(Individual population[], int i, double a0, double a1, double beta, double eps);
void simulate_ppc_beforestartday(Individual* pop, int start_day, double a0, double a1, double beta, double eps);
void simulate_ppc_afterstartday(Individual* pop, int start_day);