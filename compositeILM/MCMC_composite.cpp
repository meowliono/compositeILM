//
// Created by yirao zhang on 2023-08-09.
//
#include "MCMC_composite.h"
void MH_composite(map<int, vector<Individual>>clusters, int t_end, double samples[][4]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    samples[0][3] = 0;
    double sd_a0 = 0.01;
    double sd_a1 = 0.01;
    double sd_beta = 0.5;
    double sd_eps = 0.1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_eps);
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
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_eps) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_eps) * M_SQRT1_2));
        double nxt_prob = composite_post_dens(clusters, t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3]);
        double crt_prob = composite_post_dens(clusters, t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
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
void MH_composite_K(map<int, vector<Individual>>clusters, int t_end, double samples[][6]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    double sd_a0 = 0.01;
    double sd_a1 = 0.01;
    double sd_beta = 0.1;
    double sd_epsk = 0.001;
    for(int i=3; i<K+3; i++){
        samples[0][i] = 0;
    }
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);//generate a0
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);//generate a1
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);//generate beta
        normal_distribution<double> dis3(samples[i-1][3],sd_epsk);//generate e_k
        double nxt_sample[3+K];
        do{
            nxt_sample[0] = dis0(gen);
        }while(nxt_sample[0]<0);
        do{
            nxt_sample[1] = dis1(gen);
        }while(nxt_sample[1]<0);
        do{
            nxt_sample[2] = dis2(gen);
        }while(nxt_sample[2]<0);
        for(int i=3; i<K+3; i++){
            do{
                nxt_sample[i] = dis3(gen);
            }while(nxt_sample[i]<0);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = 1;
        double prop_nxt = 1;
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2));
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2));
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2));
        for(int j =3; j<K+3; j++){
            prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][j])/sd_epsk) * M_SQRT1_2));
            prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[j])/sd_epsk) * M_SQRT1_2));
        }
        double nxt_prob = composite_post_dens(clusters, t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], &nxt_sample[3],true);
        double crt_prob = composite_post_dens(clusters, t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], &samples[i-1][3], true);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout <<"e "<<e<< " ratio "<<ratio<<" next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            for(int j = 3; j<K+3; j++){
                samples[i][j] = nxt_sample[j];
            }
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            for(int j = 3; j<K+3; j++){
                samples[i][j] = samples[i-1][j];
            }
        }
    }
}
void MH_composite_tv(map<int, vector<Individual>>clusters, int t_end, double samples[][5]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    samples[0][3] = 0;
    samples[0][4] = 0;
    double sd_a0 = 0.008;
    double sd_a1 = 0.008;
    double sd_beta = 0.5;
    double sd_eps = 0.1;
    double sd_delta = 0.2;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_eps);
        normal_distribution<double> dis4(samples[i-1][4],sd_beta);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen), dis4(gen)};
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
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_eps) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/sd_delta) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_eps) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[4])/sd_delta) * M_SQRT1_2));
        double nxt_prob = composite_post_dens_tv(clusters, t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4]);
        double crt_prob = composite_post_dens_tv(clusters, t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<","<<samples[i][4]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<","<<samples[i][4]<<endl;
        }
    }
}
void MH_composite_bc(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][4]){//between cluster composite ILM
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.5;
    samples[0][3] = 0.001;
    double sd_a0 = 0.022;
    double sd_a1 = 0.022;
    double sd_beta = 0.2;
    double sd_eps = 0.01;
    //samples[0][4] = 1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_eps);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen)};
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
        /*
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        */
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_eps) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_eps) * M_SQRT1_2));
        double nxt_prob = composite_post_dens_bc(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3]);
        double crt_prob = composite_post_dens_bc(clusters, centroids,t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            //samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            //samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}
void MH_composite_bc_tilde(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][5]){//between cluster composite ILM
    samples[0][0] = 0.01;
    samples[0][1] = 0.01;
    samples[0][2] = 0.5;
    samples[0][3] = 0.01;
    samples[0][4] = 5;
    double sd_a0 = 0.005;
    double sd_a1 = 0.005;
    double sd_beta = 0.005;
    double sd_eps = 0.05;
    double sd_tilde = 0.05;
    //samples[0][4] = 1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_eps);
        normal_distribution<double> dis4(samples[i-1][4],sd_tilde);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen), dis4(gen)};
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
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_eps) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/sd_tilde) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_eps) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/sd_tilde) * M_SQRT1_2));
        double nxt_prob = composite_post_dens_bc_tilde(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4]);
        double crt_prob = composite_post_dens_bc_tilde(clusters, centroids, t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            //samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            //samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}
void MH_composite_bc_tilde_woeps(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][4]){//between cluster composite ILM
    samples[0][0] = 0.01;
    samples[0][1] = 0.01;
    samples[0][2] = 1;
    samples[0][3] = 1.5;
    //samples[0][4] = 1;
    double sd_a0 = 0.02;
    double sd_a1 = 0.02;
    double sd_beta = 0.4;
    double sd_tilde = 0.4;
    for(int i=1; i<20000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_tilde);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen)};
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
        /*
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        */
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_tilde) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_tilde) * M_SQRT1_2));
        double nxt_prob = composite_post_dens_bc_tilde_woeps(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3]);
        double crt_prob = composite_post_dens_bc_tilde_woeps(clusters, centroids,t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        //nxt_prob*=pow(10, 400);
        //crt_prob*=pow(10, 400);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            //samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            //samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}
void MH_composite_bc_tilde_woeps_infcentroids(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][4]){//between cluster composite ILM
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.5;
    samples[0][3] = 1.5;
    //samples[0][4] = 1;
    double sd_a0 = 0.025;
    double sd_a1 = 0.03;
    double sd_beta = 0.5;
    double sd_tilde = 0.5;
    for(int i=1; i<20000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_tilde);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen)};
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
        /*
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        */
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_tilde) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_tilde) * M_SQRT1_2));
        long double nxt_prob = composite_post_dens_bc_tilde_woeps_infcentroids(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3]);
        long double crt_prob = composite_post_dens_bc_tilde_woeps_infcentroids(clusters, centroids,t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            //samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            //samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}
void MH_composite_bc_tilde_woeps_infcentroids_delta(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][5]){//between cluster composite ILM
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.5;
    samples[0][3] = 1.5;
    samples[0][4] = 1.5;
    //samples[0][4] = 1;
    double sd_a0 = 0.025;
    double sd_a1 = 0.03;
    double sd_beta = 0.2;
    double sd_tilde = 0.2;
    double sd_delta = 0.1;
    for(int i=1; i<20000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_tilde);
        normal_distribution<double> dis4(samples[i-1][4],sd_delta);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen), dis4(gen)};
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
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_tilde) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/sd_delta) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_tilde) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[4])/sd_delta) * M_SQRT1_2));
        long double nxt_prob = composite_post_dens_bc_tilde_woeps_infcentroids_delta(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3],nxt_sample[4]);
        long double crt_prob = composite_post_dens_bc_tilde_woeps_infcentroids_delta(clusters, centroids,t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<","<<samples[i][4]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<","<<samples[i][4]<<endl;
        }
    }
}
void MH_composite_bc_tilde_woeps_infcentroids_inflated(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double samples[][4]){//between cluster composite ILM
    samples[0][0] = 0.01;
    samples[0][1] = 0.01;
    samples[0][2] = 0.5;
    samples[0][3] = 0.5;
    //samples[0][4] = 1;
    double sd_a0 = 0.02;
    double sd_a1 = 0.02;
    double sd_beta = 0.5;
    double sd_tilde = 1;
    for(int i=1; i<20000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],sd_a0);
        normal_distribution<double> dis1(samples[i-1][1],sd_a1);
        normal_distribution<double> dis2(samples[i-1][2],sd_beta);
        normal_distribution<double> dis3(samples[i-1][3],sd_tilde);
        //normal_distribution<double> dis4(samples[i-1][4],0.3);
        double nxt_sample[5] = {dis0(gen), dis1(gen), dis2(gen), dis3(gen)};
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
        /*
        while(nxt_sample[4]<0){
            nxt_sample[4] = dis4(gen);
        }
        */
        double e = (double)(rand()%100001)/100000;
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/sd_tilde) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/sd_a0) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/sd_a1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/sd_beta) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/sd_tilde) * M_SQRT1_2));
        double nxt_prob = composite_post_dens_bc_tilde_woeps_infcentroids_inflated(clusters, centroids,t_end, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3]);
        double crt_prob = composite_post_dens_bc_tilde_woeps_infcentroids_inflated(clusters, centroids,t_end, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3]);
        double ratio = prop_cur*nxt_prob/(prop_nxt*crt_prob);
        cout << "next prob " <<nxt_prob<< " current prob " <<crt_prob<< endl;
        if(ratio>=1 || (ratio<1 && e<ratio)){
            samples[i][0] = nxt_sample[0];
            samples[i][1] = nxt_sample[1];
            samples[i][2] = nxt_sample[2];
            samples[i][3] = nxt_sample[3];
            //samples[i][4] = nxt_sample[4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
        else{
            samples[i][0] = samples[i-1][0];
            samples[i][1] = samples[i-1][1];
            samples[i][2] = samples[i-1][2];
            samples[i][3] = samples[i-1][3];
            //samples[i][4] = samples[i-1][4];
            cout << "e " <<e<< "ratio " <<ratio<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<","<<samples[i][3]<<endl;
        }
    }
}