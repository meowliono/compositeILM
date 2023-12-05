#include "posterior_predictive_checks.h"
using namespace std;
int ppc_crt_time = 0;
int ppc_inf_time_tab[100][300] = {0};//a table storing each individual's infection time
int ppc_rem_time_tab[100][300] = {0};//a table storing each individual's removal time
double prob_sus_ppc(Individual population[], int id, double a0, double a1, double beta, double eps){//id - the # of the individual
    double s_i = population[id].get_sus(a0,a1);
    double sum_t_k = 0;
    double prob = 0;
    for(int i = 0; i<total_pop; i++){
        if(i == id){
            continue;
        }
        else if (i<id && population[i].cur_status == 1){
            sum_t_k += population[i].get_inf(population[i].coef[2])*population[id].kernel(population[i], beta);
        }
        else if(i>id && population[i].next_status == 1){
            sum_t_k += population[i].get_inf(population[i].coef[2])*population[id].kernel(population[i], beta);
        }
    }
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark(eps)));
    return prob;
}
double prob_sus_ppc_afterstartday(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double s_eps){//id - the # of the individual
    double s_i = s_a1*(population[id].get_sus(0,1))+s_a0;
    double sum_t_k = 0;
    double prob = 0;
    for(int j = 0; j<total_pop; j++){
        if(j == id){
            continue;
        }
        if ((ppc_inf_time_tab[i][j]!=-1) && (ppc_inf_time_tab[i][j]<=ppc_crt_time) &&(ppc_crt_time<ppc_inf_time_tab[i][j]+t_rem)){
            sum_t_k += population[j].get_inf(population[j].coef[2])*population[id].kernel(population[j], s_beta);
        }
    }
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark(s_eps)));
    return prob;
}

void status_change_ppc_beforestartday(Individual population[], int id, double a0, double a1, double beta, double eps){
    population[id].cur_status = population[id].next_status;//update current status
    if(population[id].cur_status == 0){
        double prob = prob_sus_ppc(population, id, a0, a1, beta, eps);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            population[id].next_status = 1;
            population[id].inf_time = ppc_crt_time + 1;
            for(int i = 0; i<100; i++){
                ppc_inf_time_tab[i][id] = population[id].inf_time;
            }
        }
        else{
            population[id].next_status = 0;
        }
    }
    if(population[id].cur_status == 1){
        if(ppc_crt_time+1 == population[id].inf_time + t_rem){
            population[id].next_status = 2;
            population[id].rem_time = ppc_crt_time + 1;
            for(int i = 0; i<100; i++){
                ppc_rem_time_tab[i][id] = population[id].rem_time;
            }
        }
        else if(ppc_crt_time+1 < population[id].inf_time + t_rem){
            population[id].next_status = 1;//stay infectious
        }
    }
    if(population[id].cur_status == 2){
        population[id].next_status = 2;//pre_status=2
    }
}
void status_change_ppc_afterstartday(Individual population[], int i, int id, double s_a0, double s_a1, double s_beta, double s_eps){// i is the ith row of ppc
    if(ppc_inf_time_tab[i][id]==-1){
        double prob = prob_sus_ppc_afterstartday(population, i, id, s_a0, s_a1, s_beta, s_eps);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            ppc_inf_time_tab[i][id] = ppc_crt_time + 1;
            //cout << id << "will be infected on day " << ppc_crt_time+1<< " with probability "<<prob << endl;
        }
    }
    else if(ppc_inf_time_tab[i][id]+t_rem==ppc_crt_time+1){
        ppc_rem_time_tab[i][id] = ppc_crt_time+1;

    }
}
void pop_status_change_ppc_beforestartday(Individual population[], double a0, double a1, double beta, double eps){
    for(int j = 0; j<total_pop; j++){
        status_change_ppc_beforestartday(population, j, a0, a1, beta, eps);
    }
}
void pop_status_change_ppc_afterstartday(Individual population[], int i, double s_a0, double s_a1, double s_beta, double s_eps){
    for(int j=0; j<total_pop; j++){
        status_change_ppc_afterstartday(population, i, j, s_a0, s_a1, s_beta, s_eps);
    }
}
void simulate_ppc_beforestartday(Individual* pop, int start_day, double a0, double a1, double beta, double eps) {
    for (;ppc_crt_time < start_day; ppc_crt_time++) {
        pop_status_change_ppc_beforestartday(pop, a0, a1, beta, eps);
    }
}
void simulate_ppc_afterstartday(Individual* pop, int start_day){
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
    for (; ppc_crt_time<31; ppc_crt_time++){
        for(int i = 0; i<101; i++){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(1999, 20000);
            int random_sample_id = dis(gen);
            pop_status_change_ppc_afterstartday(pop,i, samples[random_sample_id][0], samples[random_sample_id][1],samples[random_sample_id][2] , 0);
            //pop_status_change_ppc_afterstartday(pop,i, 0.03539343, 0.07298992, 1.45775779, 0);
        }
    }
}