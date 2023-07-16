#include <iostream>
using namespace std;
#include "ILM.h"
#include "generate_individuals.h"
#include "simulate.h"
#include "Inference.h"
#include "MCMC.h"
#include <string>
#include <fstream>
int main() {
    //Initialize population
    Individual pop[total_pop];
    generate(pop);
    assign_int_inf(pop,5);
    //Simulate
    simulate(pop);
    //Export infection time&removal time tables
    ofstream rec;
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/records.csv",ios::out|ios::trunc);
    rec << "id" <<","<<"time_of_infected"<<","<<"time_of_removed"<<","<<"position_x"<<","<<"position_y"<<","<<"susceptibility level"<<","<<"infectivity level"<<endl;
    for(int i = 0; i<total_pop; i++){
        rec << i <<","<<inf_time_tab[i]<<","<<rem_time_tab[i]<<","<<pop[i].position_x<<","<<pop[i].position_y<<","<<pop[i].get_sus(pop[i].coef[0], pop[i].coef[1])<<","<<pop[i].get_inf(pop[i].coef[2])<<endl;
    }
    rec.close();
    //Import csv file(epidemic data)

    //
    //cout << "----------------------"<<endl;
    /*
    double parameters[6] = {0,1,0,1,0,1};
    set_prior_Gaussian(parameters);
    RWMH(pop, samples);
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/samples.csv",ios::out|ios::trunc);
    rec << "a0" <<","<<"a1"<<","<<"n1"<<endl;
    for(int i = 0; i<10000; i++){
        rec << samples[i][0]<<","<<samples[i][1]<<","<<samples[i][2]<<endl;
    }
     */
    double parameters[6] = {0,1,0,1,0,1};
    set_prior_Gaussian(parameters);
    double samples_a0[10000];
    RWMH_a0(pop, samples_a0);
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/samples_a0.csv",ios::out|ios::trunc);
    rec << "a0" <<endl;
    for(int i = 0; i<10000; i++){
        rec << samples_a0[i] <<endl;
    }
    rec.close();
}
