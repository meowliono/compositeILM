//
// Created by yirao zhang on 2023-06-25.
//
#include "ILM.h"
using namespace std;
double a0 = 0.5;//reduce alpha
double a1 = 0.5;
double n1 = 0.5;
int beta = 1;//power law coefficient
int t_rem = 3;//infection period is 3 days
int total_pop = 100;//100 people
int inf_time_tab[100] = {0};//a table storing each individual's infection time
int rem_time_tab[100] = {0};//a table storing each individual's removal time
Individual::Individual(){
    cur_status = 0; //0: S; 1: I; 2:R
    next_status = 0;//cache previous day status
    inf_time = -1;// individual's infect time
    rem_time = -1;// individual's removal time
}
void Individual::set_sus_cov(const double& susceptibility){//susceptibility: from low to high, between 0-1
    sus_cov[0] = a0;
    sus_cov[1] = susceptibility;
}
void Individual::set_inf_cov(const double& infectivity){//infectivity: from low to high: between 0-1
    inf_cov[0] = 0;
    inf_cov[1] = infectivity;
}
void Individual::set_position(double x, double y){
    position_x = x;
    position_y = y;
}
double Individual::get_sus(){//a public interface to get individual's susceptibility
    s = sus_cov[0] + a1* sus_cov[1];
    return s;
}
double Individual::get_inf(){//a public interface to get individual's infectivity
    i = inf_cov[0]+n1*inf_cov[1];
    return i;
}
double Individual::kernel(Individual& p){//kernel function: power-law of distance
    double d = sqrt(pow((position_x - p.position_x),2) + pow((position_y - p.position_y),2));
    //distance between two individuals
    double k = pow(d,-beta);
    return k;
}
double Individual::spark(){
    double spark = 0;//default spark = 0;
    return spark;
}
double prob_sus(Individual population[], int id){//id - the # of the individual
    double s_i = population[id].get_sus();
    double sum_t_k = 0;
    double prob = 0;
    for(int i = 0; i<total_pop; i++){
        if(i == id){
            continue;
        }
        else if (i<id && population[i].cur_status == 1){
            sum_t_k += population[i].get_inf()*population[id].kernel(population[i]);
        }
        else if(i>id && population[i].next_status == 1){
            sum_t_k += population[i].get_inf()*population[id].kernel(population[i]);
        }
    }
    prob = 1-exp(-(s_i*sum_t_k+population[id].spark()));
    return prob;
}
void status_change(Individual population[], int id){
    population[id].cur_status = population[id].next_status;//update current status
    if(population[id].cur_status == 0){
        double prob = prob_sus(population, id);
        //srand((unsigned)time(NULL));
        double e = (double)(rand()%101)/100;
        if(e<=prob && prob!= 0){
            population[id].next_status = 1;
            population[id].inf_time = current_time + 1;
            inf_time_tab[id] = population[id].inf_time;
            cout << "prob=" << prob << endl;
        }
        else{
            population[id].next_status = 0;
            //cout << prob <<"no" <<e << endl;
        }
    }
    if(population[id].cur_status == 1){
        if(current_time+1 == population[id].inf_time + t_rem){
            population[id].next_status = 2;
            population[id].rem_time = current_time + 1;
            rem_time_tab[id] = population[id].rem_time;
        }
        else if(current_time+1 < population[id].inf_time + t_rem){
            population[id].next_status = 1;//stay infectious
        }
    }
    if(population[id].cur_status == 2){
        population[id].next_status = 2;//pre_status=2
    }
}
void pop_status_change(Individual population[]){
    for(int i = 0; i<total_pop; i++){
        status_change(population, i);
    }
}

