//
// Created by yirao zhang on 2023-08-09.
//
#include "MCMC_composite.h"
void MH_composite(map<int, vector<Individual>>clusters, double samples[][4]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    samples[0][3] = 0;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.01);
        normal_distribution<double> dis1(samples[i-1][1],0.01);
        normal_distribution<double> dis2(samples[i-1][2],0.5);
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
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.01) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.01) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/0.1) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.01) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.01) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/0.1) * M_SQRT1_2));
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
void MH_composite_K(map<int, vector<Individual>>clusters, double samples[][6]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    for(int i=3; i<K+3; i++){
        samples[0][i] = 0;
    }
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.01);//generate a0
        normal_distribution<double> dis1(samples[i-1][1],0.01);//generate a1
        normal_distribution<double> dis2(samples[i-1][2],0.1);//generate beta
        normal_distribution<double> dis3(samples[i-1][3],0.001);//generate e_k
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
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][0])/0.01) * M_SQRT1_2));
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][1])/0.01) * M_SQRT1_2));
        prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][2])/0.1) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[0])/0.01) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[1])/0.01) * M_SQRT1_2));
        prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[2])/0.1) * M_SQRT1_2));
        for(int j =3; j<K+3; j++){
            prop_cur *= (1-0.5 * erfc(-((0-samples[i-1][j])/0.001) * M_SQRT1_2));
            prop_nxt *= (1-0.5 * erfc(-((0-nxt_sample[j])/0.001) * M_SQRT1_2));
        }
        double ratio = prop_cur*composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], &nxt_sample[3],true)/(prop_nxt*composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], &samples[i-1][3], true));
        cout <<"e "<<e<< " ratio "<<ratio<<" next prob " <<composite_post_dens(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], &nxt_sample[3], true)<< " current prob " <<composite_post_dens(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], &samples[i-1][3], true)<< endl;
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
void MH_composite_tv(map<int, vector<Individual>>clusters, double samples[][5]){
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    samples[0][3] = 0;
    samples[0][4] = 0;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.008);
        normal_distribution<double> dis1(samples[i-1][1],0.008);
        normal_distribution<double> dis2(samples[i-1][2],0.5);
        normal_distribution<double> dis3(samples[i-1][3],0.1);
        normal_distribution<double> dis4(samples[i-1][4],0.2);
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
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.008) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.008) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/0.1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/0.2) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.008) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.008) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.5) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/0.1) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[4])/0.2) * M_SQRT1_2));
        double ratio = prop_cur*composite_post_dens_tv(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4])/(prop_nxt*composite_post_dens_tv(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4]));
        cout << "next prob " <<composite_post_dens_tv(clusters, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4])<< " current prob " <<composite_post_dens_tv(clusters, samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4])<< endl;
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
void MH_composite_bc(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, double samples[][5]){//between cluster composite ILM
    samples[0][0] = 0.001;
    samples[0][1] = 0.001;
    samples[0][2] = 0.001;
    samples[0][3] = 0.1;
    samples[0][4] = 1;
    for(int i=1; i<50000; i++){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        normal_distribution<double> dis0(samples[i-1][0],0.012);
        normal_distribution<double> dis1(samples[i-1][1],0.012);
        normal_distribution<double> dis2(samples[i-1][2],0.5);
        normal_distribution<double> dis3(samples[i-1][3],0.01);
        normal_distribution<double> dis4(samples[i-1][4],0.3);
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
        double prop_cur = (1-0.5 * erfc(-((0-samples[i-1][0])/0.012) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][1])/0.012) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][2])/0.05) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][3])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-samples[i-1][4])/0.3) * M_SQRT1_2));
        double prop_nxt = (1-0.5 * erfc(-((0-nxt_sample[0])/0.012) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[1])/0.012) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[2])/0.05) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[3])/0.02) * M_SQRT1_2))*(1-0.5 * erfc(-((0-nxt_sample[4])/0.3) * M_SQRT1_2));
        double ratio = prop_cur*composite_post_dens_bc(clusters, centroids,nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4])/(prop_nxt*composite_post_dens_bc(clusters, centroids,samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4]));
        cout << "next prob " <<composite_post_dens_bc(clusters,centroids, nxt_sample[0], nxt_sample[1], nxt_sample[2], nxt_sample[3], nxt_sample[4])<< " current prob " <<composite_post_dens_bc(clusters, centroids,samples[i-1][0], samples[i-1][1],samples[i-1][2], samples[i-1][3], samples[i-1][4])<< endl;
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