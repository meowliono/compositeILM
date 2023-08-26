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
#include <sstream>
int main() {
    //Initialize population
    Individual pop[total_pop];
    map<int, vector<Individual>> clusters;//clusters
    vector<DataPoint> centroids;
    /*
    generate_clustered(pop, 3);
    assign_int_inf(pop,5);
    //Simulate
    simulate(pop);
    //K-means clustering

    //Export infection time&removal time tables
    ofstream rec;
    rec.open("/Users/yiraozhang/CLionProjects/compositeILM/records.csv",ios::out|ios::trunc);
    rec << "id" <<","<<"time_of_infected"<<","<<"time_of_removed"<<","<<"position_x"<<","<<"position_y"<<","<<"susceptibility level"<<","<<"infectivity level"<<","<<"cluster_id"<<endl;
    for(int i = 0; i<total_pop; i++){
        rec << i <<","<<inf_time_tab[i]<<","<<rem_time_tab[i]<<","<<pop[i].position_x<<","<<pop[i].position_y<<","<<pop[i].get_sus(0,1)<<","<<pop[i].get_inf(1)<<","<<pop[i].cluster_id<<endl;
    }
    rec.close();
    */
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

        clusters = PerformClustering(pop, total_pop, K);
        centroids = calculateClusterCenters(clusters);
        ofstream rec;
        double parameters[10] = {1, 1, 1, 1, 1, 1, 1, 1,1,1};
        set_prior_parameters(parameters, 10);
        MH_composite_tv(clusters, samples);
        rec.open("/Users/yiraozhang/CLionProjects/compositeILM/samples_composite_tv.csv", ios::out | ios::trunc);
        rec << "a0" << "," << "a1" << "," << "beta"<<","<<"epsilon"<<","<<"delta"<<endl;
        for (int i = 0; i < 50000; i++) {
            rec << samples[i][0] << "," << samples[i][1] << "," << samples[i][2]<<","<<samples[i][3]<<","<<samples[i][4]<<endl;
        }
        rec.close();
}