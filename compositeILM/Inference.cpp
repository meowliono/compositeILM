//
// Created by yirao zhang on 2023-07-12.
//
#include "Inference.h"
//#include "ILM.h"
using namespace std;
/*
 * set up the prior distribution
 * default prior distribution is Gaussian
 *
 */
double log_lhd = 0;
double prior_param[6]={0};
void set_prior_Gaussian(double param[]){//mean and sd
    for(int i = 0; i<6; i++){
        prior_param[i]=param[i];
    }
}
/*
 * calculate the log likelihood
 */
double calculate_log_likelihood(Individual pop[], double a0, double a1, double beta){
    double log_lhd = 0;
    for(int t=1; t<31; t++){//time day1-30
        for(int id=0; id<total_pop; id++){
            if(pop[id].inf_time<t && pop[id].inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = pop[id].get_sus(a0, a1);
            double sum_t_k = 0;
            double prob;
            //calculate the log likelihood of infection
            for(int i=0; i<total_pop; i++){
                if(i==id){
                    continue;
                }
                if(pop[i].inf_time<=(t-1) && (t-1)<pop[i].rem_time){
                    sum_t_k += pop[i].get_inf(1)*pop[id].kernel(pop[i], beta);
                }
            }
            prob = 1-exp(-(s_i*sum_t_k+pop[id].spark()));
            if(prob!=0){
                if(pop[id].inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(pop[id].inf_time>t||pop[id].inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double post_dist_dens(Individual pop[],double a0, double a1, double beta){
    double log_lhd =0;
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[1]))*exp(-pow(((a0-prior_param[0])/prior_param[1]),2)/2));
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[3]))*exp(-pow(((a1-prior_param[2])/prior_param[3]),2)/2));
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[5]))*exp(-pow(((beta-prior_param[4])/prior_param[5]),2)/2));
    log_lhd+= calculate_log_likelihood(pop, a0, a1, beta);
    return exp(log_lhd);
}
/*
 * call the MCMC.cpp to sample from posterior distribution
 */
