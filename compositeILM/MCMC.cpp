//
// Created by yirao zhang on 2023-06-25.
//
#include "MCMC.h"
#include <armadillo>
void MH(Individual pop[], double samples[][3]){
    samples[0][0] = 0.1;
    samples[0][1] = 0.1;
    samples[0][2] = 0.1;
    for(int i=1; i<50000; i++){
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
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.5) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.5) * M_SQRT1_2));
        double ratio = prop_cur*post_dist_dens(pop, nxt_sample[0], nxt_sample[1], nxt_sample[2])/(prop_nxt*post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2]));
        cout << "next prob " <<post_dist_dens(pop, nxt_sample[0], nxt_sample[1], nxt_sample[2])<< " current prob " <<post_dist_dens(pop, samples[i-1][0], samples[i-1][1],samples[i-1][2])<< endl;
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
            cout << "e " <<e<< "ratio " <<ratio<< " "<<nxt_sample[0] <<" "<<nxt_sample[1]<<" "<<nxt_sample[2]<<endl;
        }
    }
}
void MH_composite(map<int, vector<Individual>>clusters, double samples[][4]){
    samples[0][0] = 0.1;
    samples[0][1] = 0.1;
    samples[0][2] = 0.1;
    samples[0][3] = 0.1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.5);
        normal_distribution<double> dis1(samples[i-1][1],0.5);
        normal_distribution<double> dis2(samples[i-1][2],0.5);
        normal_distribution<double> dis3(samples[i-1][3],0.5);
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
        double ratio = composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3])/composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        cout << "next prob " <<composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3])<< " current prob " <<composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3])<< endl;
        if(ratio>=1 && e<1){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<" "<<samples[i][3]<<endl;
        }
        else if(ratio<1 && e<ratio){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<" "<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<" "<<samples[i][3]<<endl;
        }
    }
}

void leapfrog(Individual pop[], double z[], double r[], double e, int L){//z-parameters(a0, a1, beta), r-anxillary variable, e-step size, L-number of iterations
    //double E = -log(post_dist_dens(pop, z[0], z[1], z[3])); //momentum
    double znxt[3];
    double rnxt[3];
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<double> dist(0,1);
    double rand = dist(gen);
    if(rand > 0.5){
        e = -e;
    }
}

void HMC(Individual pop[], double samples[][3]){

}

void MH_log(Individual pop[], double samples[][3]){
    /*
    samples[0][0] = 1e-10;
    samples[0][1] = 1e-10;
    samples[0][2] = 1e-10;
    */
    samples[0][0] = 0.01;
    samples[0][1] = 0.01;
    samples[0][2] = 0.01;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(log(samples[i-1][0]),0.1);
        normal_distribution<double> dis1(log(samples[i-1][1]),0.1);
        normal_distribution<double> dis2(log(samples[i-1][2]),0.1);
        double nxt_sample[3] = {dis0(gen), dis1(gen), dis2(gen)};
        double e = (double)(rand()%100001)/100000;
        double ratio = post_dist_dens_logprior(pop, nxt_sample[0], nxt_sample[1], nxt_sample[2])/post_dist_dens_logprior(pop, log(samples[i-1][0]), log(samples[i-1][1]),log(samples[i-1][2]));
        cout << post_dist_dens_logprior(pop,nxt_sample[0], nxt_sample[1], nxt_sample[2])<<endl;
        cout << ratio <<endl;
        if(ratio>=1 && e<1){
            samples[i][0] = exp(nxt_sample[0]);
            samples[i][1] = exp(nxt_sample[1]);
            samples[i][2] = exp(nxt_sample[2]);
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<endl;
        }
        else if(ratio<1 && e<ratio){
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