//
// Created by yirao zhang on 2023-06-25.
//
#include "ILM.h"
double prob_sus(Individual population[], int id){//id - the # of the individual
    double s_i;
    double sum_t_k = 0;
    double prob = 0;
    int total = sizeof(population)/sizeof(population[0]);
    for(int i = 0; i<total; i++){
        if(i == id){
            continue;
        }
        if(population[i].pre_status == 1){
            sum_t_k += population[i].get_inf()*population[id].kernel(population[i]);
        }
    }
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark()));
}
void status_change(Individual population[], int id){
    if(population[id].status == 0){
        double prob = prob_sus(population, id);
        std::bernoulli_distribution b(prob);
        std::default_random_engine e;
        if(b(e)){
            population[id].pre_status = population[id].status; //pre_status=0
            population[id].status = 1;
            population[id].inf_time = current_time;
        }
    }
    if(population[id].status == 1){
        if(current_time == population[id].inf_time + t_rem){
            population[id].pre_status = population[id].status;//pre_status=1
            population[id].status = 2;
            population[id].rem_time = current_time;
        }
        else if(current_time < population[id].inf_time + t_rem){
            population[id].pre_status = population[id].status;//pre_status=1
        }
    }
    if(population[id].status == 2){
        population[id].pre_status = population[id].status;//pre_status=2
    }
}
void pop_status_change(Individual population[]){
    int total = sizeof(population)/sizeof(population[0]);
    for(int i = 0; i<total; i++){
        status_change(population, i);
    }
}

