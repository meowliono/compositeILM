#include <iostream>
using namespace std;
#include "ILM.h"
#include "generate_individuals.h"
#include "simulate.h"
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
    rec << "id" <<","<<"time_of_infected"<<","<<"time_of_removed"<<","<<"position_x"<<","<<"position_y"<<endl;
    for(int i = 0; i<100; i++){
        rec << i <<","<<inf_time_tab[i]<<","<<rem_time_tab[i]<<","<<pop[i].position_x<<","<<pop[i].position_y<<endl;
    }
    rec.close();
}
