//
// Created by yirao zhang on 2023-06-25.
//
#include "MCMC.h"
#include <armadillo>
using namespace std;
using namespace arma;
void MH(Individual pop[], double samples[][3]){
    /*
    samples[0][0] = 1e-10;
    samples[0][1] = 1e-10;
    samples[0][2] = 1e-10;
    */
    samples[0][0] = 0.1;
    samples[0][1] = 0.1;
    samples[0][2] = 0.1;
    mat cov = { {2, -0.1, 0.1},
              {-0.1, 2, 0.1},
                {0.1,0.1,2}};
    for(int i=1; i<50000; i++){
        /*
         * block-wise update
         * generating log(a0), log(a1), log(beta)
        vec cur_sample_log = {log(samples[i-1][0]), log(samples[i-1][1]), log(samples[i-1][2])};
        vec nxt_sample_log = mvnrnd(cur_sample_log, cov);
        double e = (double)(rand()%100001)/100000;
        double ratio = post_dist_dens(pop, exp(nxt_sample_log(0)), exp(nxt_sample_log(1)), exp(nxt_sample_log(2)))/post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2]);
        cout << "next prob " <<post_dist_dens(pop, exp(nxt_sample_log(0)), exp(nxt_sample_log(1)), exp(nxt_sample_log(2)))<< " current prob " <<post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2])<< endl;
        if(ratio>=1 && e<1){
            samples[i][0] = exp(nxt_sample_log(0));
            samples[i][1] = exp(nxt_sample_log(1));
            samples[i][2] = exp(nxt_sample_log(2));
            cout << "e " <<e<< "ratio " <<ratio<< " "<<exp(nxt_sample_log(0)) <<" "<<exp(nxt_sample_log(1))<<" "<<exp(nxt_sample_log(2))<<endl;
        }
        else if(ratio<1 && e<ratio){
            samples[i][0] = exp(nxt_sample_log(0));
            samples[i][1] = exp(nxt_sample_log(1));
            samples[i][2] = exp(nxt_sample_log(2));
            cout << "e " <<e<< "ratio " <<ratio<< " "<<exp(nxt_sample_log(0)) <<" "<<exp(nxt_sample_log(1))<<" "<<exp(nxt_sample_log(2))<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<exp(nxt_sample_log(0)) <<" "<<exp(nxt_sample_log(1))<<" "<<exp(nxt_sample_log(2))<<endl;
        }
        */
        /* independently update a0, a1, beta
        */
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.5);
        normal_distribution<double> dis1(samples[i-1][1],0.5);
        normal_distribution<double> dis2(samples[i-1][2],0.5);
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
        double ratio = post_dist_dens(pop, nxt_sample[0], nxt_sample[1], nxt_sample[2])/post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2]);
        cout << "next prob " <<post_dist_dens(pop, nxt_sample[0], nxt_sample[1], nxt_sample[2])<< " current prob " <<post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2])<< endl;
        if(ratio>=1 && e<1){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample[0] <<" "<<nxt_sample[1]<<" "<<nxt_sample[2]<<endl;
        }
        else if(ratio<1 && e<ratio){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample[0] <<" "<<nxt_sample[1]<<" "<<nxt_sample[2]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample[0] <<" "<<nxt_sample[1]<<" "<<nxt_sample[2]<<endl;

        }
    }
}
void MH_a0(Individual pop[], double samples_a0[]){
    for(int i=1; i<10000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples_a0[i-1],0.8);
        double nxt_sample = dis0(gen);
        while(nxt_sample<0){
            nxt_sample = dis0(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double ratio = post_dist_dens(pop, nxt_sample, 0.5, 0.5)/post_dist_dens(pop, samples_a0[i-1], 0.5, 0.5);
        cout << "next prob " <<post_dist_dens(pop, nxt_sample, 0.5, 0.5)<< " current prob " <<post_dist_dens(pop, samples_a0[i-1], 0.5, 0.5)<< endl;
        if(ratio>=1 && e<1){
            samples_a0[i] = nxt_sample;
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample <<endl;
        }
        else if(ratio<1 && e<ratio){
            samples_a0[i] = nxt_sample;
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample <<endl;
        }
        else{
            samples_a0[i] = samples_a0[i-1];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples_a0[i-1] <<endl;
        }
    }
}