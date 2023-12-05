#include "posterior_predictive_checks.h"
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
            //cout << "individual"<<id<<"will be infected on day"<<(current_time+1)<<"the prob is"<<prob<<endl;
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
void status_change_ppc_afterstartday(Individual population[], int i, int id, double a0, double a1, double beta, double eps){// i is the ith row of ppc
    population[id].cur_status = population[id].next_status;//update current status
    if(population[id].cur_status == 0){
        double prob = prob_sus_ppc(population, id, a0, a1, beta, eps);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            population[id].next_status = 1;
            population[id].inf_time = ppc_crt_time + 1;
            ppc_inf_time_tab[i][id] = population[id].inf_time;
            //cout << "individual"<<id<<"will be infected on day"<<(current_time+1)<<"the prob is"<<prob<<endl;
        }
        else{
            population[id].next_status = 0;
        }
    }
    if(population[id].cur_status == 1){
        if(ppc_crt_time+1 == population[id].inf_time + t_rem){
            population[id].next_status = 2;
            population[id].rem_time = ppc_crt_time + 1;
            ppc_rem_time_tab[i][id] = population[id].rem_time;
        }
        else if(ppc_crt_time+1 < population[id].inf_time + t_rem){
            population[id].next_status = 1;//stay infectious
        }
    }
    if(population[id].cur_status == 2){
        population[id].next_status = 2;//pre_status=2
    }
}
void pop_status_change_ppc_beforestartday(Individual population[], double a0, double a1, double beta, double eps){
    for(int i = 0; i<total_pop; i++){
        status_change_ppc_beforestartday(population, i, a0, a1, beta, eps);
    }
}
void pop_status_change_ppc_afterstartday(Individual population[], int i, double a0, double a1, double beta, double eps){
    for(int j=0; j<total_pop; j++){
        status_change_ppc_afterstartday(population, i, j, a0, a1, beta, eps);
    }
}
void simulate_ppc_beforestartday(Individual* pop, int start_day, double a0, double a1, double beta, double eps) {
    for (int i = 0; i < start_day; i++) {
        pop_status_change_ppc_beforestartday(pop, a0, a1, beta, eps);
        ppc_crt_time++;
    }
}
void simulate_ppc_afterstartday(Individual* pop, int start_day){
    for (int i = start_day; i<31; i++){
        for(int j = 0; j<100; j++){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(999, 10000);
            int random_sample_id = dis(gen);
            pop_status_change_ppc_afterstartday(pop,j,samples[random_sample_id][0] , samples[random_sample_id][1], samples[random_sample_id][2], 0);
        }
        ppc_crt_time++;
    }
}