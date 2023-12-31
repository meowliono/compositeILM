//
// Created by yirao zhang on 2023-07-12.
//
#include "Inference.h"
using namespace std;
/*
 * set up the prior distribution
 * default prior distribution is Gaussian
 *
 */
vector<double> prior_param(8,0);
void set_prior_parameters(double param[], int num){
    for(int i = 0; i<num; i++){
        prior_param[i]=param[i];
    }
}
/*
 * calculate the log likelihood
 */
double calculate_log_likelihood(Individual pop[], int t_end, double a0, double a1, double beta){
    double log_lhd = 0;
    for(int t=1; t<=t_end; t++){//time day1-30
        for(int id=0; id<total_pop; id++){
            if((pop[id].inf_time<t) && (pop[id].inf_time!=-1)){//if the individual is already infected, skip
                continue;
            }
            double s_i = pop[id].get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(int i=0; i<total_pop; i++){
                if(i==id){
                    continue;
                }
                if((pop[i].inf_time!=-1) && (pop[i].inf_time<=(t-1)) && ((t-1)<pop[i].rem_time||pop[i].rem_time==-1)){
                    sum_t_k += pop[i].get_inf(1)*pop[id].kernel(pop[i], beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k));
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
double calculate_log_likelihood_logprior(Individual pop[], int t_end, double ta0, double ta1, double tbeta){//ta0-log(a0); ta1-log(a1); ta2-log(a2)
    double log_lhd = 0;
    for(int t=1; t<t_end; t++){//time day1-30
        for(int id=0; id<total_pop; id++){
            if(pop[id].inf_time<t && pop[id].inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = pop[id].get_sus(exp(ta0), exp(ta1));
            double sum_t_k = 0;
            double prob;
            //calculate the log likelihood of infection
            for(int i=0; i<total_pop; i++){
                if(i==id){
                    continue;
                }
                if(pop[i].inf_time<=(t-1) && ((t-1)<pop[i].rem_time||pop[i].rem_time==-1)){
                    sum_t_k += pop[i].get_inf(1)*pop[id].kernel(pop[i], exp(tbeta));
                }
            }
            prob = 1-exp(-(s_i*sum_t_k+pop[id].spark(pop[id].coef[4])));
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
double post_dist_dens(Individual pop[],int t_end, double a0, double a1, double beta){
    double log_lhd =0;
    /* normal distribution as prior
     *
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[1]))*exp(-pow(((a0-prior_param[0])/prior_param[1]),2)/2))+log(2);
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[3]))*exp(-pow(((a1-prior_param[2])/prior_param[3]),2)/2))+log(2);
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[5]))*exp(-pow(((beta-prior_param[4])/prior_param[5]),2)/2))+log(2);
    log_lhd+= calculate_log_likelihood(pop, a0, a1, beta);
     */
    /*
     * gamma distribution as prior
     */
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    log_lhd+= calculate_log_likelihood(pop, t_end, a0, a1, beta);
    return exp(log_lhd);
}
double post_dist_dens_logprior(Individual pop[],int t_end, double ta0, double ta1, double tbeta){
    double log_lhd =0;
    /* normal distribution as prior
     *
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[1]))*exp(-pow(((a0-prior_param[0])/prior_param[1]),2)/2))+log(2);
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[3]))*exp(-pow(((a1-prior_param[2])/prior_param[3]),2)/2))+log(2);
    log_lhd+= log((1.0/(sqrt(2*M_PI)*prior_param[5]))*exp(-pow(((beta-prior_param[4])/prior_param[5]),2)/2))+log(2);
    log_lhd+= calculate_log_likelihood(pop, a0, a1, beta);
     */
    /*
     * gamma distribution as prior
     */
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, exp(ta0)))+ta0;
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, exp(ta1)))+ta1;
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, exp(tbeta)))+tbeta;
    log_lhd += calculate_log_likelihood_logprior(pop, t_end, ta0, ta1, tbeta);
    return exp(log_lhd);
}
/*
 * calculate the composite likelihood
 */
