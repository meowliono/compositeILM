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
double log_lhd = 0;
double prior_param[6]={0};
void set_prior_parameters(double param[]){
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
            //calculate the log likelihood of infection
            for(int i=0; i<total_pop; i++){
                if(i==id){
                    continue;
                }
                if(pop[i].inf_time<=(t-1) && (t-1)<pop[i].rem_time){
                    sum_t_k += pop[i].get_inf(1)*pop[id].kernel(pop[i], beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+pop[id].spark(pop[id].coef[4])));
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
double calculate_log_likelihood_logprior(Individual pop[], double ta0, double ta1, double tbeta){//ta0-log(a0); ta1-log(a1); ta2-log(a2)
    double log_lhd = 0;
    for(int t=1; t<31; t++){//time day1-30
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
                if(pop[i].inf_time<=(t-1) && (t-1)<pop[i].rem_time){
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
double post_dist_dens(Individual pop[],double a0, double a1, double beta){
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
    log_lhd+= calculate_log_likelihood(pop, a0, a1, beta);
    return exp(log_lhd);
}
double post_dist_dens_logprior(Individual pop[],double ta0, double ta1, double tbeta){
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
    log_lhd += log(boost::math::pdf(gammaDist1, exp(ta0)))+2*ta0;
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, exp(ta1)))+2*ta1;
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, exp(tbeta)))+2*tbeta;
    log_lhd += calculate_log_likelihood_logprior(pop, ta0, ta1, tbeta);
    return exp(log_lhd);
}
/*
 * calculate the composite likelihood
 */
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps){
    double log_lhd = 0;
    for(int t=1; t<31; t++){//time day1-30
        for(Individual& person : clusters[K]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            double prob;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[K]){
                if(&other==&person){
                    continue;
                }
                if(other.inf_time<=(t-1) && (t-1)<other.rem_time){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            prob = 1-exp(-(s_i*sum_t_k+person.spark(eps)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double composite_log_likelihood(map<int, vector<Individual>> clusters, double a0, double a1, double beta, double eps){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood(clusters, i, a0, a1, beta, eps);
    }
    return log_lhd;
}
double composite_post_dens(map<int, vector<Individual>>clusters, double a0, double a1, double beta, double eps){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, eps));
    log_lhd+= composite_log_likelihood(clusters, a0, a1, beta, eps);
    return exp(log_lhd);
}