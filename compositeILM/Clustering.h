
#ifndef COMPOSITEILM_CLUSTERING_H
#define COMPOSITEILM_CLUSTERING_H
#include "ILM.h"
#include <map>
#include <vector>
using namespace std;
#endif //COMPOSITEILM_CLUSTERING_H
void clustering_PAM(Individual pop[], int K);
void clustering_Kmeans(Individual pop[], int K);
map<int, vector<Individual>>PerformClustering(Individual pop[], int total_pop, int K);