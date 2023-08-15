//
// Created by yirao zhang on 2023-08-09.
//
#include "MCMC_composite.h"
void MH_composite(map<int, vector<Individual>>clusters, double samples[][4]){
    samples[0][0] = 0.1;
    samples[0][1] = 0.1;
    samples[0][2] = 0.1;
    samples[0][3] = 0.1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.2);
        normal_distribution<double> dis1(samples[i-1][1],0.2);
        normal_distribution<double> dis2(samples[i-1][2],1);
        normal_distribution<double> dis3(samples[i-1][3],0.1);
        double nxt_sample[4] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen)};
        while(nxt_sample[0]<0){
            nxt_sample[0] = dis0(gen);
        }
        while(nxt_sample[1]<0){
            nxt_sample[1] = dis1(gen);
        }
        while(nxt_sample[2]<0){
            nxt_sample[2] = dis2(gen);
        }
        while(nxt_sample[3]<0){
            nxt_sample[3] = dis3(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/0.5) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/0.5) * M_SQRT1_2));
        double ratio = prop_cur*composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3])/(prop_nxt*composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]));
        cout << "next prob " <<composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3])<< " current prob " <<composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3])<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}