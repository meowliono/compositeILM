//
// Created by yirao zhang on 2023-06-26.
#include "generate_individuals.h"
#include <armadillo>
using namespace arma;

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
        //pop[i].set_inf_cov((double)(rand()%101)/100);
        pop[i].set_inf_cov(1);
    }
}
void generate_clustered(Individual pop[], int K){//generate the clustered population
    //armadillo library to generate bivariate Gaussian distribution
    for(int i=0; i<(K-1); i++){
        float x = static_cast<float>((i+1))/(K+1);
        vec2 mean = {5*x,5*x};
        mat22 cov = {{0.2, 0}, {0, 0.2}};
        arma_rng::set_seed_random();
        // Generate a random sample from the bivariate Gaussian distribution
        for(int j=0; j<(total_pop/K); j++){
            vec2 sample = mvnrnd(mean, cov);
            inf_time_tab[i*(total_pop/K)+j] = -1;
            rem_time_tab[i*(total_pop/K)+j] = -1;
            pop[i*(total_pop/K)+j].set_position(sample(0), sample(1));
            pop[i*(total_pop/K)+j].set_sus_cov((double)(rand()%101)/100);
            //pop[i*(total_pop/K)+j].set_inf_cov((double)(rand()%101)/100);
            pop[i*(total_pop/K)+j].set_inf_cov(1);
            pop[i*(total_pop/K)+j].cluster_id = i;
        }
    }
    float x = static_cast<float>(K)/(K+1);
    vec2 mean = {5*x,5*x};
    mat22 cov = {{0.1, 0}, {0, 0.1}};
    arma_rng::set_seed_random();
    for(int j=(K-1)*(total_pop/K); j<total_pop; j++){
        vec2 sample = mvnrnd(mean, cov);
        inf_time_tab[j] = -1;
        rem_time_tab[j] = -1;
        pop[j].set_position(sample(0), sample(1));
        pop[j].set_sus_cov((double)(rand()%101)/100);
        //pop[j].set_inf_cov((double)(rand()%101)/100);
        pop[j].set_inf_cov(1);
        pop[j].cluster_id = K-1;
    }
}
void assign_int_inf(Individual pop[], int id){//identify the 'patient 0'
    pop[id].cur_status = 1;
    pop[id].next_status = 1;
    pop[id].inf_time = 0;
    inf_time_tab[id] = 0;
}
