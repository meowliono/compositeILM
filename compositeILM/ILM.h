//
// Created by yirao zhang on 2023-06-25.
//

#ifndef COMPOSITEILM_ILM_H
#define COMPOSITEILM_ILM_H

#endif //COMPOSITEILM_ILM_H
#include <math.h>
#include <random>
int a0 = 1;
int a1 = 1;
int n1 = 1;
int beta = 2;//power law coefficient
int t_rem = 3;//infection period is 3 days
int current_time = 0;//current time
class Individual{//define individual class
private:
    double s;// individual's susceptibility
    double sus_cov[2];//only two covariates: a0-intercept; a1-susceptibility
    double i;// individual's infectivity
    double inf_cov[2];//only two covariates: a0-intercept; a1-infectivity
    double position_x;// individual's position x axis
    double position_y;// individual's position y axis
public:
    int status = 0; //0: S; 1: I; 2:R
    int pre_status = 0;//cache previous day status
    int inf_time = 0;// individual's infect time
    int rem_time = 0;// individual's removal time
public:
    /*
     set_sus_cov: set susceptibility covariates
     get_sus: get susceptibility-s
     set_inf_cov: set infectivity covariates
     get_inf: get infectivity-i
     */
    void set_sus_cov(const double& susceptibility){//susceptibility: from low to high: 0, 1, 2
        sus_cov[0] = a0;
        sus_cov[1] = susceptibility;
    }
    void set_inf_cov(const double& infectivity){//infectivity: from low to high: 0, 1, 2
        inf_cov[0] = 0;
        inf_cov[1] = infectivity;
    }
    void set_position(double x, double y){
        position_x = x;
        position_y = y;
    }
    double get_sus(){//a public interface to get individual's susceptibility
        s = sus_cov[0] + a1* sus_cov[1];
        return s;
    }
    double get_inf(){//a public interface to get individual's infectivity
        i = inf_cov[0]+n1*inf_cov[1];
        return i;
    }
    double kernel(Individual& p){//kernel function: power-law of distance
        double d = sqrt(pow((position_x - p.position_x),2) + pow((position_y - p.position_y),2));
        //distance between two individuals
        double k = pow(d,-beta);
        return k;
    }
    double spark(){
        double spark = 0;//default spark = 0;
        return spark;
    }
};