//
// Created by yirao zhang on 2023-06-26.
#include "generate_individuals.h"
#include <armadillo>
using namespace std;

/* randomly generate 100 individuals in a 100*100 area
 * randomly assign their susceptibility and infectivity
 * randomly decide which individuals are infectious at time 0
*/
void generate(Individual pop[]){//the total number of individuals
    srand((unsigned)time(NULL));
    for(int i = 0; i<total_pop; i++){
        inf_time_tab[i] = -1;
        rem_time_tab[i] = -1;
        for(int r = 0; r<100; r++){
            ppc_inf_time_tab[r][i] = -1;
            ppc_rem_time_tab[r][i] = -1;
        }
        pop[i].set_position(rand()%11+((double)(rand()%101)/100), rand()%11+((double)(rand()%101)/100));
        /*
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        gamma_distribution<double> dis(2, 1);
        pop[i].set_sus_cov(dis(gen)/10);
         */
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        gamma_distribution<double> dis1(2, 1);
        gamma_distribution<double> dis2(3, 1);
        double x = dis1(gen);
        double y = dis2(gen);
        pop[i].set_sus_cov(x/(x+y));
        pop[i].set_inf_cov(1);
    }
}
void generate_clustered(Individual pop[], int K){//generate the clustered population
    //armadillo library to generate bivariate Gaussian distribution
    for(int i=0; i<(K-1); i++){
        arma::vec2 mean = {(double)((rand()%3001)/100),(double)((rand()%3001)/100)};
        arma::mat22 cov = {{8, 0}, {0, 8}};
        arma::arma_rng::set_seed_random();
        // Generate a random sample from the bivariate Gaussian distribution
        for(int j=0; j<(total_pop/K); j++){
            arma::vec2 sample = mvnrnd(mean, cov);
            inf_time_tab[i*(total_pop/K)+j] = -1;
            rem_time_tab[i*(total_pop/K)+j] = -1;
            for(int r = 0; r<100; r++){
                ppc_inf_time_tab[r][i*(total_pop/K)+j] = -1;
                ppc_rem_time_tab[r][i*(total_pop/K)+j] = -1;
            }
            pop[i*(total_pop/K)+j].set_position(sample(0), sample(1));
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            default_random_engine gen(seed);
            gamma_distribution<double> dis1(2, 1);
            gamma_distribution<double> dis2(3, 1);
            double x = dis1(gen);
            double y = dis2(gen);
            pop[i*(total_pop/K)+j].set_sus_cov(x/(x+y));
            //pop[i*(total_pop/K)+j].set_sus_cov((double)(rand()%101)/100);
            pop[i*(total_pop/K)+j].set_inf_cov(1);
            pop[i*(total_pop/K)+j].cluster_id = i;
        }
    }
    arma::vec2 mean = {(double)((rand()%3001)/100),(double)((rand()%3001)/100)};
    arma::mat22 cov = {{8, 0}, {0, 8}};
    arma::arma_rng::set_seed_random();
    for(int j=(K-1)*(total_pop/K); j<total_pop; j++){
        arma::vec2 sample = mvnrnd(mean, cov);
        inf_time_tab[j] = -1;
        rem_time_tab[j] = -1;
        for(int r = 0; r<100; r++){
            ppc_inf_time_tab[r][j] = -1;
            ppc_rem_time_tab[r][j] = -1;
        }
        pop[j].set_position(sample(0), sample(1));
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        gamma_distribution<double> dis1(2, 1);
        gamma_distribution<double> dis2(3, 1);
        double x = dis1(gen);
        double y = dis2(gen);
        pop[j].set_sus_cov(x/(x+y));
        //pop[j].set_sus_cov((double)(rand()%101)/100);
        pop[j].set_inf_cov(1);
        pop[j].cluster_id = K-1;
    }
}
void generate_clustered_random(Individual pop[], int K){//generate the clustered population
    //armadillo library to generate bivariate Gaussian distribution
    for(int i=0; i<(K-1); i++){
        arma::vec2 mean = {(double)((rand()%3001)/100),(double)((rand()%3001)/100)};
        arma::mat22 cov = {{8, 0}, {0, 8}};
        arma::arma_rng::set_seed_random();
        // Generate a random sample from the bivariate Gaussian distribution
        for(int j=0; j<(total_pop/K); j++){
            arma::vec2 sample = mvnrnd(mean, cov);
            inf_time_tab[i*(total_pop/K)+j] = -1;
            rem_time_tab[i*(total_pop/K)+j] = -1;
            for(int r = 0; r<100; r++){
                ppc_inf_time_tab[r][i*(total_pop/K)+j] = -1;
                ppc_rem_time_tab[r][i*(total_pop/K)+j] = -1;
            }
            pop[i*(total_pop/K)+j].set_position(sample(0), sample(1));
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            default_random_engine gen(seed);
            gamma_distribution<double> dis1(2, 1);
            gamma_distribution<double> dis2(3, 1);
            double x = dis1(gen);
            double y = dis2(gen);
            pop[i*(total_pop/K)+j].set_sus_cov(x/(x+y));
            //pop[i*(total_pop/K)+j].set_sus_cov((double)(rand()%101)/100);
            pop[i*(total_pop/K)+j].set_inf_cov(1);
            //pop[i*(total_pop/K)+j].cluster_id = i;
        }
    }
    arma::vec2 mean = {(double)((rand()%3001)/100),(double)((rand()%3001)/100)};
    arma::mat22 cov = {{8, 0}, {0, 8}};
    arma::arma_rng::set_seed_random();
    for(int j=(K-1)*(total_pop/K); j<total_pop; j++){
        arma::vec2 sample = mvnrnd(mean, cov);
        inf_time_tab[j] = -1;
        rem_time_tab[j] = -1;
        for(int r = 0; r<100; r++){
            ppc_inf_time_tab[r][j] = -1;
            ppc_rem_time_tab[r][j] = -1;
        }
        pop[j].set_position(sample(0), sample(1));
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);
        gamma_distribution<double> dis1(2, 1);
        gamma_distribution<double> dis2(3, 1);
        double x = dis1(gen);
        double y = dis2(gen);
        pop[j].set_sus_cov(x/(x+y));
        //pop[j].set_sus_cov((double)(rand()%101)/100);
        pop[j].set_inf_cov(1);
        //pop[j].cluster_id = K-1;
    }
}
void assign_int_inf(Individual pop[], int id){//identify the 'patient 0'
    pop[id].cur_status = 1;
    pop[id].next_status = 1;
    pop[id].inf_time = 0;
    inf_time_tab[id] = 0;
    //ppc_inf_time_tab[id] = 0;
    for(int r = 0; r<100; r++){
        ppc_inf_time_tab[r][id] = 0;
    }
}
