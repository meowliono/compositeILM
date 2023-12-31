//
// Created by yirao zhang on 2023-06-25.
//
#include "MCMC.h"
void MH(Individual pop[], int t_end, double samples[][3]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    for(int i=1; i<20000; i++){
        /* independently update a0, a1, beta
         */
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.02);
        normal_distribution<double> dis1(samples[i-1][1],0.02);
        normal_distribution<double> dis2(samples[i-1][2],0.3);
        double nxt_sample[3] = {dis0(gen), dis1(gen), dis2(gen)};
        while(nxt_sample[0]<0){
            nxt_sample[0] = dis0(gen);
        }
        while(nxt_sample[1]<0){
            nxt_sample[1] = dis1(gen);
        }
        while(nxt_sample[2]<0){
            nxt_sample[2] = dis2(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.3) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.3) * M_SQRT1_2));
        double cur_prob = post_dist_dens(pop, t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2]);
        double nxt_prob = post_dist_dens(pop, t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*cur_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<cur_prob << endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample[0] <<" "<<nxt_sample[1]<<" "<<nxt_sample[2]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i-1][0] <<" "<<samples[i-1][1]<<" "<<samples[i-1][2]<<endl;
        }
    }
}

void MH_log(Individual pop[], int t_end, double samples[][3]){
    /*
    samples[0][0] = 1e-10;
    samples[0][1] = 1e-10;
    samples[0][2] = 1e-10;
    */
    samples[0][0] = 0.1;
    samples[0][1] = 0.1;
    samples[0][2] = 0.1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(log(samples[i-1][0]),0.4);
        normal_distribution<double> dis1(log(samples[i-1][1]),0.4);
        normal_distribution<double> dis2(log(samples[i-1][2]),0.4);
        double nxt_sample[3] = {dis0(gen), dis1(gen), dis2(gen)};
        double e = (double)(rand()%100001)/100000;
        double ratio = post_dist_dens_logprior(pop, t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2])/post_dist_dens_logprior(pop, t_end, log(samples[i-1][0]), log(samples[i-1][1]),log(samples[i-1][2]));
        if(ratio>=1 ||(ratio<1 && e<ratio)){
            samples[i][0] = exp(nxt_sample[0]);
            samples[i][1] = exp(nxt_sample[1]);
            samples[i][2] = exp(nxt_sample[2]);
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<endl;
        }
    }
}