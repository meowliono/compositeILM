/*
 * clustering methods including PAM, K-means...
 */
#include "Clustering.h"
using namespace std;
void clustering_PAM(Individual pop[], int K){

}

void clustering_Kmeans(Individual pop[], int k){
    double centroids[k][2];
    //initialize centroids
    int seed = (int)time(NULL);
    srand(seed);
    for(int i = 0; i<k; i++){
        int c = rand()%total_pop;
        pop[c].cluster_id = i;
        centroids[i][0] = pop[c].position_x;
        centroids[i][1] = pop[c].position_y;
    }
    bool stop = false;
    int count;
    while(!stop || count >1000){
        stop = true;
        //assign classes
        for(int i = 0; i<total_pop; i++){
            double min_dis = sqrt(pow((centroids[0][0] - pop[i].position_x),2) + pow((centroids[0][1] - pop[i].position_y), 2));
            for(int j= 1; j<k; j++){
                //calculate the distance between each individual and k centroids
                double dis = sqrt(pow((centroids[i][0] - pop[i].position_x),2) + pow((centroids[i][1] - pop[i].position_y), 2));
                if(dis<min_dis){
                    min_dis = dis;
                    if(pop[i].cluster_id == j){
                        stop = false;
                    }
                    pop[i].cluster_id = j;
                }
            }
        }
        //calculate new centroids
        if (stop == false){
            int mark_clear[k];//mark whether the cluster has been visited, if not, clear the current centroids
            for(int i = 0; i<total_pop; i++){
                if(!mark_clear[pop[i].cluster_id]){
                    centroids[pop[i].cluster_id][0]=pop[i].position_x;
                    centroids[pop[i].cluster_id][1]=pop[i].position_y;
                    mark_clear[pop[i].cluster_id] = 1;
                }
                else{
                    centroids[pop[i].cluster_id][0]+=pop[i].position_x;
                    centroids[pop[i].cluster_id][1]+=pop[i].position_y;
                }
            }
        }
        count++;
    }
}
//K-means algorithm in opencv

//assign population to k centroids
map<int, vector<Individual>>PerformClustering(Individual pop[], int total_pop, int K){
    map<int, vector<Individual>> clusters;
    for(int i = 0; i<total_pop; i++){
        clusters[pop[i].cluster_id].push_back(pop[i]);
    }
    return clusters;
}
