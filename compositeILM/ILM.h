//
// Created by yirao zhang on 2023-06-25.
//
#ifndef COMPOSITEILM_ILM_H
#define COMPOSITEILM_ILM_H
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
extern float beta;//power law coefficient
extern int t_rem;//infection period is 3 days
extern int current_time;//current time
extern int total_pop;//total population
extern int inf_time_tab[];//a table store individuals' infection time
extern int rem_time_tab[];//a table store individuals' removal time
extern double samples[][3];
class Individual{//define individual class
private:
    double s;// individual's susceptibility
    double sus_cov[2];//only two covariates: a0-intercept; a1-susceptibility
    double i;// individual's infectivity
    double inf_cov[2];//only two covariates: a0-intercept; a1-infectivity
public:
    double coef[3];//a0,a1,n1
    int cur_status; //0: S; 1: I; 2:R store current day status
    int next_status;//store the calculated next day status
    int inf_time;// individual's infection time
    int rem_time;// individual's removal time
    double position_x;// individual's position x axis
    double position_y;// individual's position y axis
public:
    /*
     set_sus_cov: set susceptibility covariates
     get_sus: get susceptibility-s
     set_inf_cov: set infectivity covariates
     get_inf: get infectivity-i
     */
    Individual();
    void set_sus_cov(const double& susceptibility);
    void set_inf_cov(const double& infectivity);
    void set_position(double x, double y);
    double get_sus(double a0, double a1);
    double get_inf(double n1);
    double kernel(Individual& p);
    double spark();
};
double prob_sus(Individual population[], int id);
void status_change(Individual population[], int id);
void pop_status_change(Individual population[]);
#endif //COMPOSITEILM_ILM_H
