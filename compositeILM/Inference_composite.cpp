//
// Created by yirao zhang on 2023-08-09.
//
#include "Inference_composite.h"
using namespace std;
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int K, double a0, double a1, double beta, double eps){
    double log_lhd = 0;
    for(int t=1; t<31; t++){//time day1-30
        for(Individual& person : clusters[K]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[K]){
                if(&other==&person){
                    continue;
                }
                if(other.inf_time<=(t-1) && (t-1)<other.rem_time){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(eps)));
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