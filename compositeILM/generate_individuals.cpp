//
// Created by yirao zhang on 2023-06-26.
#include "generate_individuals.h"
/* randomly generate 100 individuals in a 100*100 area
 * randomly assign their susceptibility and infectivity
 * randomly decide which individuals are infectious at time 0
*/
void generate(Individual pop[]){//the total number of individuals
    srand((unsigned)time(NULL));
    for(int i = 0; i<total_pop; i++){
        inf_time_tab[i] = -1;
        rem_time_tab[i] = -1;
        pop[i].set_position(rand()%101+((double)(rand()%101)/100), rand()%101+((double)(rand()%101)/100));
        pop[i].set_sus_cov((double)(rand()%101)/100);
        pop[i].set_inf_cov((double)(rand()%101)/100);
    }
}
void assign_int_inf(Individual pop[], int id){//identify the 'patient 0'
    pop[id].cur_status = 1;
    pop[id].next_status = 1;
    pop[id].inf_time = 0;
    inf_time_tab[id] = 0;
}
