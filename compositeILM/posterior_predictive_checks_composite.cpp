//
// Created by yirao zhang on 2023-10-14.
//
#include "posterior_predictive_checks_composite.h"
using namespace std;
double prob_sus_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, int id, double s_a0, double s_a1, double s_beta, double t_beta){//id - the # of the individual
    double s_i = s_a1*(population[id].get_sus(0,1))+s_a0;
    double sum_t_k = 0;
    double prob = 0;
    double sum_term = 0;//sum of distance multiplied by number of infections in each cluster
    vector<double> distances(K,0);
    vector<double> num_by_clusters(K, 0);
    for(int j = 0; j<total_pop; j++){
        if(j == id){
            continue;
        }
        if ((ppc_inf_time_tab[i][j]!=-1) && (ppc_inf_time_tab[i][j]<=ppc_crt_time)&&(ppc_crt_time<(ppc_inf_time_tab[i][j]+t_rem))){
            if(population[j].cluster_id == population[id].cluster_id){
                sum_t_k += population[j].get_inf(population[j].coef[2])*population[id].kernel(population[j], s_beta);
            }
            else{
                num_by_clusters[population[j].cluster_id] += 1;
            }
        }
    }
    for (int k = 0; k< K; k++){
        if(k!=population[id].cluster_id){
            if(num_by_clusters[k]!=0){
                distances[k] = sqrt(pow((centroids[population[id].cluster_id].x-centroids[k].x),2) + pow((centroids[population[id].cluster_id].y-centroids[k].y), 2));
                sum_term += num_by_clusters[k] * pow(distances[k], -t_beta);
            }
        }
    }
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark(s_i*sum_term)));
    return prob;
}
double prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta){//id - the # of the individual
    double s_i = s_a1*(population[id].get_sus(0,1))+s_a0;
    double sum_t_k = 0;
    double prob = 0;
    double sum_term = 0;//sum of distance multiplied by number of infections in each cluster
    vector<DataPoint> centroids(K,{0,0});
    vector<int> num_by_clusters(K, 0);
    for(int j = 0; j<total_pop; j++){
        if(j == id){
            continue;
        }
        else if ((ppc_inf_time_tab[i][j]>-1) && (ppc_inf_time_tab[i][j]<=ppc_crt_time) &&(ppc_crt_time<(ppc_inf_time_tab[i][j]+t_rem))){
            if(population[j].cluster_id == population[id].cluster_id){
                sum_t_k += population[j].get_inf(population[j].coef[2])*population[id].kernel(population[j], s_beta);
            }
            else{
                num_by_clusters[population[j].cluster_id] += 1;
                centroids[population[j].cluster_id].x += population[j].position_x;
                centroids[population[j].cluster_id].y += population[j].position_y;
            }
        }
    }
    for (int k = 0; k< K; k++){
        if(k!=population[id].cluster_id){
            if(num_by_clusters[k]!=0){
                centroids[k].x = centroids[k].x/num_by_clusters[k];
                centroids[k].y = centroids[k].y/num_by_clusters[k];
                double distance = sqrt(pow((centroids[population[id].cluster_id].x-centroids[k].x),2) + pow((centroids[population[id].cluster_id].y-centroids[k].y), 2));
                sum_term += num_by_clusters[k] * pow(distance, -t_beta);
            }
        }
    }
    //cout << "sum_t_k " << sum_t_k <<" sum_term " << sum_term << endl;
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark(s_i*sum_term)));
    return prob;
}
double prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta){//id - the # of the individual
    double s_i = s_a1*(population[id].get_sus(0,1))+s_a0;
    double sum_t_k = 0;
    double prob = 0;
    double sum_term = 0;//sum of distance multiplied by number of infections in each cluster
    vector<DataPoint> centroids(K,{0,0});
    vector<int> num_by_clusters(K, 0);
    for(int j = 0; j<total_pop; j++){
        if(j == id){
            continue;
        }
        else if ((ppc_inf_time_tab[i][j]>-1) && (ppc_inf_time_tab[i][j]<=ppc_crt_time) &&(ppc_crt_time<(ppc_inf_time_tab[i][j]+t_rem))){
            if(population[j].cluster_id == population[id].cluster_id){
                sum_t_k += population[j].get_inf(population[j].coef[2])*population[id].kernel(population[j], s_beta);
            }
            else{
                num_by_clusters[population[j].cluster_id] += 1;
                centroids[population[j].cluster_id].x += population[j].position_x;
                centroids[population[j].cluster_id].y += population[j].position_y;
            }
        }
    }
    for (int k = 0; k< K; k++){
        if(k!=population[id].cluster_id){
            if(num_by_clusters[k]!=0){
                centroids[k].x = centroids[k].x/num_by_clusters[k];
                centroids[k].y = centroids[k].y/num_by_clusters[k];
                double distance = sqrt(pow((centroids[population[id].cluster_id].x-centroids[k].x),2) + pow((centroids[population[id].cluster_id].y-centroids[k].y), 2));
                sum_term += pow(num_by_clusters[k],s_delta) * pow(distance, -t_beta);
            }
        }
    }
    //cout << "sum_t_k " << sum_t_k <<" sum_term " << sum_term << endl;
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark(s_i*sum_term)));
    return prob;
}
void status_change_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, int id, double s_a0, double s_a1, double s_beta, double t_beta){// i is the ith row of ppc
    if(ppc_inf_time_tab[i][id]==-1){
        double prob = prob_sus_ppc_afterstartday_bc_tilde_woeps(population, centroids,i, id, s_a0, s_a1, s_beta,t_beta);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            ppc_inf_time_tab[i][id] = ppc_crt_time + 1;
            //cout << "in the " << i << " th posterior" << id << " individual was infected with prob = " << prob << endl;
        }
    }
    else if(ppc_inf_time_tab[i][id]+t_rem==(ppc_crt_time+1)){
        ppc_rem_time_tab[i][id] = ppc_crt_time+1;
    }
}
void status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta){// i is the ith row of ppc
    if(ppc_inf_time_tab[i][id]==-1){
        double prob = prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids(population, i, id, s_a0, s_a1, s_beta,t_beta);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            ppc_inf_time_tab[i][id] = ppc_crt_time + 1;
            //cout << "in the " << i << " th posterior" << id << " individual was infected with prob = " << prob << endl;
        }
    }
    else if(ppc_inf_time_tab[i][id]+t_rem==ppc_crt_time+1){
        ppc_rem_time_tab[i][id] = ppc_crt_time+1;
    }
}
void status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta){// i is the ith row of ppc
    if(ppc_inf_time_tab[i][id]==-1){
        double prob = prob_sus_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(population, i, id, s_a0, s_a1, s_beta,t_beta, s_delta);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            ppc_inf_time_tab[i][id] = ppc_crt_time + 1;
            //cout << "in the " << i << " th posterior" << id << " individual was infected with prob = " << prob << endl;
        }
    }
    else if(ppc_inf_time_tab[i][id]+t_rem==ppc_crt_time+1){
        ppc_rem_time_tab[i][id] = ppc_crt_time+1;
    }
}
void pop_status_change_ppc_afterstartday_bc_tilde_woeps(Individual population[], vector<DataPoint> centroids, int i, double s_a0, double s_a1, double s_beta, double t_beta){
    for(int j=0; j<total_pop; j++){
        status_change_ppc_afterstartday_bc_tilde_woeps(population, centroids, i, j, s_a0, s_a1, s_beta, t_beta);
    }
}
void pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual population[], int i, double s_a0, double s_a1, double s_beta, double t_beta){
    for(int j=0; j<total_pop; j++){
        status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(population, i, j, s_a0, s_a1, s_beta, t_beta);
    }
}
void pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual population[], int i, double s_a0, double s_a1, double s_beta, double t_beta, double s_delta){
    for(int j=0; j<total_pop; j++){
        status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(population, i, j, s_a0, s_a1, s_beta, t_beta, s_delta);
    }
}
void simulate_ppc_afterstartday_bc_tilde_woeps(Individual* pop, vector<DataPoint> centroids, int start_day, int end_day){
    for(int t = 0;t < start_day;t++){
        for(int i = 1; i < total_pop; i++){
            if(pop[i].inf_time==t){
                for(int j = 1; j<100; j++){
                    ppc_inf_time_tab[j][i] = t;
                }
            }
        }
    }
    ppc_crt_time = start_day-1;
    for (; ppc_crt_time<end_day; ppc_crt_time++){
        for(int i = 0; i<100; i++){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(1999, 20000);
            int random_sample_id = dis(gen);
            pop_status_change_ppc_afterstartday_bc_tilde_woeps(pop, centroids, i, samples[random_sample_id][0], samples[random_sample_id][1], samples[random_sample_id][2] , samples[random_sample_id][3]-1);
            //pop_status_change_ppc_afterstartday(pop,i, 0.03539343, 0.07298992, 1.45775779, 0);
        }
    }
}
void simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids(Individual* pop, int start_day, int end_day){
    for(int t = 0;t < start_day;t++){
        for(int i = 1; i < total_pop; i++){
            if(pop[i].inf_time==t){
                for(int j = 1; j<100; j++){
                    ppc_inf_time_tab[j][i] = t;
                }
            }
        }
    }
    ppc_crt_time = start_day-1;
    for (; ppc_crt_time<end_day; ppc_crt_time++){
        for(int i = 0; i<100; i++){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(1999, 20000);
            int random_sample_id = dis(gen);
            pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids(pop,i, samples[random_sample_id][0],samples[random_sample_id][1],samples[random_sample_id][2] , samples[random_sample_id][3]);
            //pop_status_change_ppc_afterstartday(pop,i, 0.03539343, 0.07298992, 1.45775779, 0);es[random_sample_id][3]
        }
    }
}
void simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual* pop, int start_day, int end_day){
    for(int t = 0;t < start_day;t++){
        for(int i = 1; i < total_pop; i++){
            if(pop[i].inf_time==t){
                for(int j = 1; j<100; j++){
                    ppc_inf_time_tab[j][i] = t;
                }
            }
        }
    }
    ppc_crt_time = start_day-1;
    for (; ppc_crt_time<end_day; ppc_crt_time++){
        for(int i = 0; i<100; i++){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(1999, 20000);
            int random_sample_id = dis(gen);
            pop_status_change_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(pop,i, samples[random_sample_id][0],samples[random_sample_id][1],samples[random_sample_id][2] , samples[random_sample_id][3], samples[random_sample_id][4]);
            //pop_status_change_ppc_afterstartday(pop,i, 0.03539343, 0.07298992, 1.45775779, 0);es[random_sample_id][3]
        }
    }
}
void partially_simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(Individual* pop, int stepsize){
    int last_matched = 31;
    for(int t=31; t>=0; t--){
        cout << "current time is" << t << endl;
        if(t==last_matched-stepsize){
            cout << "matched is " << t << endl;
            if(stepsize == 1){
                simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(pop, t+1, last_matched);
            }
            else{
                simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(pop, t + 1, last_matched - 1);
                }
            last_matched = t;
        }
    }
}
void partially_simulate_ppc_afterstartday_bc_tilde_woeps(Individual* pop, vector<DataPoint> centroids, const vector<int>& records){
    for (; ppc_crt_time<31; ){
        auto it = find(records.begin(), records.end(), ppc_crt_time);
        if(it!=records.end()){//from it to it+2 were observed
            for(int i = 0; i < total_pop; i++){
                if(ppc_crt_time<=pop[i].inf_time && pop[i].inf_time<=ppc_crt_time+2){
                    for (int j = 1; j < 100; j++) {
                        ppc_inf_time_tab[j][i] = pop[i].inf_time;
                    }
                }
                if(ppc_crt_time>pop[i].inf_time && pop[i].inf_time!=-1){
                    for(int j = 1; j<100; j++){
                        if(ppc_inf_time_tab[j][i]<ppc_crt_time && ppc_inf_time_tab[j][i]!=-1){
                            continue;
                        }
                        else{
                            ppc_inf_time_tab[j][i]=-100;//already infected individuals can be tested and set it to a very small value
                        }
                    }
                }
                if(pop[i].inf_time==-1 || pop[i].inf_time>ppc_crt_time+2){
                    for(int j = 1; j<100; j++){
                        ppc_inf_time_tab[j][i] = -1;
                    }
                }
            }
            ppc_crt_time += 2;
        }
        else{
            for(int i = 0; i<100; i++){
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<int> dis(1999, 20000);
                int random_sample_id = dis(gen);
                pop_status_change_ppc_afterstartday_bc_tilde_woeps(pop, centroids, i, samples[random_sample_id][0],samples[random_sample_id][1],samples[random_sample_id][2] , samples[random_sample_id][3]);
                //pop_status_change_ppc_afterstartday(pop,i, 0.03539343, 0.07298992, 1.45775779, 0);es[random_sample_id][3]
            }
            ppc_crt_time++;
        }
    }
}