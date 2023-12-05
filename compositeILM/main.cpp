#include <iostream>
using namespace std;
#include "ILM.h"
#include "generate_individuals.h"
#include "simulate.h"
#include "Inference.h"
#include "MCMC.h"
#include <string>
#include <fstream>
#include "Clustering.h"
#include "HamiltonianMC.h"
#include "Inference_composite.h"
#include "MCMC_composite.h"
#include "posterior_predictive_checks.h"
#include "posterior_predictive_checks_composite.h"
#include <sstream>
int main() {
    //Initialize population
    Individual pop[total_pop];
    map<int, vector<Individual>> clusters;//clusters
    vector<DataPoint> centroids;
    ofstream rec;
    generate_clustered(pop, K);
    assign_int_inf(pop,5);

    //Simulate
    simulate(pop);
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/records.csv",ios::out|ios::trunc);
    rec << "id" <<","<<"time_of_infected"<<","<<"time_of_removed"<<","<<"position_x"<<","<<"position_y"<<","<<"susceptibility level"<<","<<"infectivity level"<<","<<"cluster_id"<<endl;
    for(int i = 0; i<total_pop; i++){
        rec << i <<","<<inf_time_tab[i]<<","<<rem_time_tab[i]<<","<<pop[i].position_x<<","<<pop[i].position_y<<","<<pop[i].get_sus(0,1)<<","<<pop[i].get_inf(1)<<","<<pop[i].cluster_id<<endl;
    }
    rec.close();

    /*
    //Import csv file(epidemic data)
    fstream epidata;
    string line;
    epidata.open("/Users/yiraozhang/CLionProjects/compositeILM/records.csv", ios::in);
    getline(epidata, line);//skip the header line
    for (int i = 0; i < total_pop; i++) {
        getline(epidata, line);
        stringstream ss(line);
        string element;
        getline(ss, element, ',');//get id
        getline(ss, element, ',');//get infection time
        pop[i].inf_time = stoi(element);
        inf_time_tab[i] = pop[i].inf_time;
        getline(ss, element, ',');//get remove time
        pop[i].rem_time = stoi(element);
        rem_time_tab[i] = pop[i].rem_time;
        getline(ss, element, ',');//get position_x
        pop[i].position_x = stod(element);
        getline(ss, element, ',');//get position_y
        pop[i].position_y = stod(element);
        getline(ss, element, ',');//get susceptibility level(between 0-1)
        pop[i].set_sus_cov(stod(element));
        getline(ss, element, ',');//get infectivity level(1 by default)
        pop[i].set_inf_cov(stod(element));
        getline(ss, element, ',');
        pop[i].cluster_id = stoi(element);
    }
    */
    /*
    //MCMC sampling
    clusters = PerformClustering(pop, total_pop, K);
    centroids = calculateClusterCenters(clusters);
    double parameters[8] = {1, 1, 1, 1, 2.5, 1, 3, 1};
    set_prior_parameters(parameters, 8);
    MH_composite_bc_tilde_woeps(clusters, centroids, 31, samples);
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/samples_composite_bc_tilde_woeps.csv", ios::out | ios::trunc);
    rec << "a0" << "," << "a1" << "," << "beta"<<","<<"tilde_beta"<<endl;
    for (int i = 0; i < 20000; i++) {
        rec << samples[i][0] << "," << samples[i][1] << "," << samples[i][2]<<","<<samples[i][3]<<endl;
    }
    rec.close();
    */
    /*
    //Posterior predictive distribution
    fstream sample;
    string line2;
    sample.open("/Users/yiraozhang/CLionProjects/compositeILM/samples_composite_bc_tilde_woeps_infcentroids_delta.csv", ios::in);
    getline(sample, line2);//skip the header line
    for (int i = 0; i < 20000; i++) {
        getline(sample, line2);
        stringstream ss(line2);
        string element;
        getline(ss, element, ',');//get a0
        samples[i][0]= stod(element);
        getline(ss, element, ',');//get a1
        samples[i][1]= stod(element);
        getline(ss, element, ',');//get beta
        samples[i][2]= stod(element);
        getline(ss, element, ',');//get t_beta
        samples[i][3]= stod(element);
        getline(ss, element, ',');//get s_delta
        samples[i][4]= stod(element);
    }
    clusters = PerformClustering(pop, total_pop, K);
    centroids = calculateClusterCenters(clusters);
    //vector<int> observed = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
    int step_size = 3;
    simulate_ppc_afterstartday_bc_tilde_woeps_infcentroids_delta(pop, 11, 31);
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/posterior_predictive_checks_inftime.csv", ios::out | ios::trunc);
    for (int i = 0; i < 100; i++) {
        for(int j=0; j<300; j++){
            rec << ppc_inf_time_tab[i][j] << ",";
        }
        rec << endl;
    }
    rec.close();
    */
}