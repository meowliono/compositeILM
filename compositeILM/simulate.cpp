//
// Created by yirao zhang on 2023-06-27.
//
#include "simulate.h"
#include "ILM.h"
using namespace std;
int current_time = 0;
void simulate(Individual* pop){
    for (int i = 0; i<30; i++){
        pop_status_change(pop);
        cout << "current time = " << current_time << endl;
        current_time++;
    }
}
