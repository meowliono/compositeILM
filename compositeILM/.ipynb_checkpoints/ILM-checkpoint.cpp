//
// Created by yirao zhang on 2023-06-25.
//
#include "ILM.h"
using namespace std;
/*double a0 = 0.5;//reduce alpha
 * double a1 = 0.5;
 * double n1 = 0.5;
 */
float beta = 2;//power law coefficient
int t_rem = 3;//infection period is 3 days
int total_pop = 100;//100 people
int inf_time_tab[100] = {0};//a table storing each individual's infection time
int rem_time_tab[100] = {0};//a table storing each individual's removal time
double samples[50000][3] = {0};
Individual::Individual(){
    cur_status = 0; //0: S; 1: I; 2:R
    next_status = 0;//cache previous day status
    inf_time = -1;// individual's infect time
    rem_time = -1;// individual's removal time
    coef[0] = 2;//a0
    coef[1] = 2;//a1
    coef[2] = 1;//n1
    coef[3] = 1.5;//beta
    cluster_id = 0;
}
void Individual::set_sus_cov(const double& susceptibility){//susceptibility: from low to high, between 0-1
    sus_cov[0] = coef[0];
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
double Individual::get_sus(double a0, double a1){//a public interface to get individual's susceptibility
    s = a0 + a1* sus_cov[1];
    return s;
}
double Individual::get_inf(double n1){//a public interface to get individual's infectivity
    i = inf_cov[0]+n1*inf_cov[1];
    return i;
}
double Individual::kernel(Individual& p, double beta){//kernel function: power-law of distance
    double d = sqrt(pow((position_x - p.position_x),2) + pow((position_y - p.position_y), 2));
    //distance between two individuals
    double k = pow(d,-beta);
    return k;
}
double Individual::spark(){
    double spark = 0;//default spark = 0;
    return spark;
}
double prob_sus(Individual population[], int id){//id - the # of the individual
    double s_i = population[id].get_sus(population[id].coef[0],population[id].coef[1]);
    double sum_t_k = 0;
    double prob = 0;
    for(int i = 0; i<total_pop; i++){
        if(i == id){
            continue;
        }
        else if (i<id && population[i].cur_status == 1){
            sum_t_k += population[i].get_inf(population[i].coef[2])*population[id].kernel(population[i], population[id].coef[3]);
        }
        else if(i>id && population[i].next_status == 1){
            sum_t_k += population[i].get_inf(population[i].coef[2])*population[id].kernel(population[i], population[id].coef[3]);
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
        double e = (double)(rand()%100001)/100000;
        if(e<=prob && prob!= 0){
            population[id].next_status = 1;
            population[id].inf_time = current_time + 1;
            inf_time_tab[id] = population[id].inf_time;
            //cout << "individual"<<id<<"will be infected on day"<<(current_time+1)<<"the prob is"<<prob<<endl;
        }
        else{
            population[id].next_status = 0;
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

