//
// Created by yirao zhang on 2023-08-09.
//
#include "Inference_composite.h"
using namespace std;
double cluster_log_likelihood(map<int, vector<Individual>> clusters, int t_end, int id, double a0, double a1, double beta, double eps){
    double log_lhd = 0;
    for(int t=1; t<t_end; t++){//time day1-t_end
        for(Individual& person : clusters[id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection within cluster id
            for(Individual& other : clusters[id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1)) && ((t-1)<other.rem_time||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(eps)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_tv(map<int, vector<Individual>> clusters, int t_end, int id, double a0, double a1, double beta, double eps, double delta){
    double log_lhd = 0;
    for(int t=1; t<t_end; t++){//time day1-t_end
        int num = 0;//record the number of infections
        for(int i=0; i<total_pop; i++){
            if((inf_time_tab[i]!=-1)&&(inf_time_tab[i]<=(t-1)) &&(t-1<rem_time_tab[i]||rem_time_tab[-1])){
                num++;
            }
        }
        for(Individual& person : clusters[id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1)) && ((t-1)<other.rem_time)||other.rem_time==-1){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(eps*pow(num,delta))));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_bc(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end,int id, double a0, double a1, double beta, double eps){
    double log_lhd = 0;
    //calculate the distances between clusters
    vector<double> distances(K,0);
    for(int i=0; i<K; i++){
        if(i!=id){
            distances[i] = sqrt(pow((centroids[id].x-centroids[i].x),2) + pow((centroids[id].y-centroids[i].y), 2));
        }
        else{
            distances[i] = 0;
        }
    }
    for(int t=1; t<t_end; t++){//time day1-30
        vector<int> num_by_clusters(K,0);//record the number of infected in each cluster
        double sum_term;
        for(int i=0; i<K; i++){
            for(Individual& person: clusters[i]){
                if((person.inf_time!=-1)&&(person.inf_time<=(t-1)) && ((t-1)<person.rem_time||person.rem_time==-1)){
                    num_by_clusters[i]++;
                }
            }
            if(i!=id){
                sum_term += num_by_clusters[i] * pow(distances[i], -beta);
            }
        }
        for(Individual& person : clusters[id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1)) && ((t-1)<other.rem_time||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(eps*sum_term)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_bc_tilde(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, int id, double a0, double a1, double beta, double eps, double tilde_beta){
    double log_lhd = 0;
    //calculate the distances between clusters
    vector<double> distances(K,0);
    for(int i=0; i<K; i++){
        if(i!=id){
            distances[i] = sqrt(pow((centroids[id].x-centroids[i].x),2) + pow((centroids[id].y-centroids[i].y), 2));
        }
        else{
            distances[i] = 0;
        }
    }
    for(int t=1; t<t_end; t++){//time day1-30
        vector<int> num_by_clusters(K,0);//record the number of infected in each cluster
        double sum_term;
        for(int i=0; i<K; i++){
            for(Individual& person: clusters[i]){
                if((person.inf_time!=-1)&&(person.inf_time<=(t-1)) && ((t-1)<person.rem_time||person.rem_time==-1)){
                    num_by_clusters[i]++;
                }
            }
            if(i!=id){
                sum_term += num_by_clusters[i] * pow(distances[i], -tilde_beta);
            }
        }
        for(Individual& person : clusters[id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1)) && ((t-1)<other.rem_time||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(eps*sum_term)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_bc_tilde_woeps(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, int id, double a0, double a1, double beta, double tilde_beta){
    double log_lhd = 0;
    //calculate the distances between clusters
    vector<double> distances(K,0);
    for(int i=0; i<K; i++){
        if(i!=id){
            distances[i] = sqrt(pow((centroids[id].x-centroids[i].x),2) + pow((centroids[id].y-centroids[i].y), 2));
        }
        else{
            distances[i] = 0;
        }
    }
    for(int t=1; t<t_end; t++){//time day1-30
        vector<int> num_by_clusters(K,0);//record the number of infected in each cluster
        double sum_term;
        for(int i=0; i<K; i++){
            for(Individual& person: clusters[i]){
                if((person.inf_time!=-1)&&(person.inf_time<=(t-1)) && ((t-1)<person.rem_time||person.rem_time==-1)){
                    num_by_clusters[i]++;
                }
            }
            if(i!=id){
                sum_term += num_by_clusters[i] * pow(distances[i], -tilde_beta);
            }
        }
        for(Individual& person : clusters[id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1)) && ((t-1)<other.rem_time||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(s_i*sum_term)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_bc_tilde_woeps_infcentroids(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, int cluster_id, double a0, double a1, double beta, double tilde_beta){
    double log_lhd = 0;
    //calculate the distances between clusters
    vector<double> distances(K,0);
    for(int t=1; t<t_end; t++){//time day1-30
        vector<int> num_by_clusters(K,0);//record the number of infected in each cluster
        double sum_term;
        for(int i=0; i<K; i++){
            if(i == cluster_id){
                continue;
            }
            double x;//centroids of infected individuals in cluster i
            double y;//centroids of infected individuals in cluster i
            for(Individual& person: clusters[i]){
                if((person.inf_time!=-1)&&(person.inf_time<=(t-1)) && ((t-1)<person.rem_time||person.rem_time==-1)){
                    num_by_clusters[i]++;
                    x+=person.position_x;
                    y+=person.position_y;
                }
            }
            if(num_by_clusters[i]!=0){
                x = x/num_by_clusters[i];
                y = y/num_by_clusters[i];
                distances[i] = sqrt(pow((centroids[cluster_id].x-x),2) + pow((centroids[cluster_id].y-y), 2));
                sum_term += num_by_clusters[i] * pow(distances[i], -tilde_beta);
            }
        }
        for(Individual& person : clusters[cluster_id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[cluster_id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1))&& (((t-1)<other.rem_time)||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(s_i*sum_term)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double cluster_log_likelihood_bc_tilde_woeps_infcentroids_delta(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, int cluster_id, double a0, double a1, double beta, double tilde_beta, double delta){
    double log_lhd = 0;
    //calculate the distances between clusters
    vector<double> distances(K,0);
    for(int t=1; t<t_end; t++){//time day1-30
        vector<int> num_by_clusters(K,0);//record the number of infected in each cluster
        double sum_term;
        for(int i=0; i<K; i++){
            if(i == cluster_id){
                continue;
            }
            double x;//centroids of infected individuals in cluster i
            double y;//centroids of infected individuals in cluster i
            for(Individual& person: clusters[i]){
                if((person.inf_time!=-1)&&(person.inf_time<=(t-1)) && ((t-1)<person.rem_time||person.rem_time==-1)){
                    num_by_clusters[i]++;
                    x+=person.position_x;
                    y+=person.position_y;
                }
            }
            if(num_by_clusters[i]!=0){
                x = x/num_by_clusters[i];
                y = y/num_by_clusters[i];
                distances[i] = sqrt(pow((centroids[cluster_id].x-x),2) + pow((centroids[cluster_id].y-y), 2));
                sum_term += pow(num_by_clusters[i],delta)* pow(distances[i], -tilde_beta);
            }
        }
        for(Individual& person : clusters[cluster_id]){
            if(person.inf_time<t && person.inf_time!=-1){//if the individual is already infected, skip
                continue;
            }
            double s_i = person.get_sus(a0, a1);
            double sum_t_k = 0;
            //calculate the log likelihood of infection
            for(Individual& other : clusters[cluster_id]){
                if(&other==&person){
                    continue;
                }
                if((other.inf_time!=-1)&&(other.inf_time<=(t-1))&& (((t-1)<other.rem_time)||other.rem_time==-1)){
                    sum_t_k += other.get_inf(1)*person.kernel(other, beta);
                }
            }
            double prob = 1-exp(-(s_i*sum_t_k+person.spark(s_i*sum_term)));
            if(prob!=0){
                if(person.inf_time==t){
                    log_lhd+=log(prob);
                    //cout << "individual"<<id<<"is infected on day"<<t<<"the prob is"<<prob<<endl;
                }
                if(person.inf_time>t||person.inf_time==-1){
                    log_lhd+=log(1-prob);
                }
            }
        }
    }
    return log_lhd;
}
double composite_log_likelihood(map<int, vector<Individual>> clusters, int t_end, double a0, double a1, double beta, double eps){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood(clusters, t_end,i, a0, a1, beta, eps);
    }
    return log_lhd;
}
double composite_log_likelihood(map<int, vector<Individual>> clusters, int t_end, double a0, double a1, double beta, double eps[],bool ClusterDepedent){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood(clusters, t_end, i, a0, a1, beta, eps[i]);
    }
    return log_lhd;
}
double composite_log_likelihood_tv(map<int, vector<Individual>> clusters, int t_end, double a0, double a1, double beta, double eps, double delta){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_tv(clusters, t_end,i, a0, a1, beta, eps, delta);
    }
    return log_lhd;
}
double composite_log_likelihood_bc(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double eps){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_bc(clusters, centroids, t_end,i, a0, a1, beta, eps);
    }
    return log_lhd;
}
double composite_log_likelihood_bc_tilde(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double eps, double tilde_beta){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_bc_tilde(clusters, centroids, t_end,i, a0, a1, beta, eps, tilde_beta);
    }
    return log_lhd;
}
double composite_log_likelihood_bc_tilde_woeps(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_bc_tilde_woeps(clusters, centroids, t_end,i, a0, a1, beta, tilde_beta);
    }
    return log_lhd;
}
double composite_log_likelihood_bc_tilde_woeps_infcentroids(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_bc_tilde_woeps_infcentroids(clusters, centroids, t_end,i, a0, a1, beta, tilde_beta);
    }
    return log_lhd;
}
double composite_log_likelihood_bc_tilde_woeps_infcentroids_delta(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta, double delta){
    double log_lhd;
    for(int i=0; i<K; i++){
        log_lhd+= cluster_log_likelihood_bc_tilde_woeps_infcentroids_delta(clusters, centroids, t_end,i, a0, a1, beta, tilde_beta, delta);
    }
    return log_lhd;
}
double composite_log_likelihood_bc_tilde_woeps_infcentroids_inflated(map<int, vector<Individual>> clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd;
    int random_i = rand()%K;
    log_lhd = K*cluster_log_likelihood_bc_tilde_woeps_infcentroids(clusters, centroids, t_end,random_i, a0, a1, beta, tilde_beta);
    return log_lhd;
}
double composite_post_dens(map<int, vector<Individual>>clusters, int t_end, double a0, double a1, double beta, double eps){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, eps));
    log_lhd+= composite_log_likelihood(clusters, t_end, a0, a1, beta, eps);
    return exp(log_lhd);
}
double composite_post_dens(map<int, vector<Individual>>clusters, int t_end, double a0, double a1, double beta, double eps[], bool ClusterDependent){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    for(int i = 0; i<K; i++){
        boost::math::gamma_distribution<> gammaDist(prior_param[6+i], prior_param[7+i]);
        log_lhd += log(boost::math::pdf(gammaDist, eps[i]));
    }
    log_lhd+= composite_log_likelihood(clusters, t_end, a0, a1, beta, eps, true);
    return exp(log_lhd);
}
double composite_post_dens_tv(map<int, vector<Individual>>clusters, int t_end, double a0, double a1, double beta, double eps, double delta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, eps));
    boost::math::gamma_distribution<> gammaDist5(prior_param[8], prior_param[9]);
    log_lhd += log(boost::math::pdf(gammaDist5, delta));
    log_lhd+= composite_log_likelihood_tv(clusters, t_end, a0, a1, beta, eps, delta);
    return exp(log_lhd);
}
double composite_post_dens_bc(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double eps){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, eps));
    log_lhd+= composite_log_likelihood_bc(clusters, centroids, t_end, a0, a1, beta, eps);
    return exp(log_lhd);
}
double composite_post_dens_bc_tilde(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double eps, double tilde_beta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, eps));
    boost::math::gamma_distribution<> gammaDist5(prior_param[8], prior_param[9]);
    log_lhd += log(boost::math::pdf(gammaDist5, eps));
    log_lhd+= composite_log_likelihood_bc_tilde(clusters, centroids, t_end, a0, a1, beta, eps, tilde_beta);
    return exp(log_lhd);
}
double composite_post_dens_bc_tilde_woeps(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, tilde_beta));
    log_lhd+= composite_log_likelihood_bc_tilde_woeps(clusters, centroids, t_end, a0, a1, beta, tilde_beta);
    return exp(log_lhd);
}
double composite_post_dens_bc_tilde_woeps_infcentroids(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, tilde_beta));
    log_lhd+= composite_log_likelihood_bc_tilde_woeps_infcentroids(clusters, centroids, t_end, a0, a1, beta, tilde_beta);
    return exp(log_lhd);
}
double composite_post_dens_bc_tilde_woeps_infcentroids_delta(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta, double delta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, tilde_beta));
    boost::math::gamma_distribution<> gammaDist5(prior_param[8], prior_param[9]);
    log_lhd += log(boost::math::pdf(gammaDist5, delta));
    log_lhd+= composite_log_likelihood_bc_tilde_woeps_infcentroids_delta(clusters, centroids, t_end, a0, a1, beta, tilde_beta, delta);
    return exp(log_lhd);
}
double composite_post_dens_bc_tilde_woeps_infcentroids_inflated(map<int, vector<Individual>>clusters, vector<DataPoint> centroids, int t_end, double a0, double a1, double beta, double tilde_beta){
    double log_lhd =0;
    boost::math::gamma_distribution<> gammaDist1(prior_param[0], prior_param[1]);
    log_lhd += log(boost::math::pdf(gammaDist1, a0));
    boost::math::gamma_distribution<> gammaDist2(prior_param[2], prior_param[3]);
    log_lhd += log(boost::math::pdf(gammaDist2, a1));
    boost::math::gamma_distribution<> gammaDist3(prior_param[4], prior_param[5]);
    log_lhd += log(boost::math::pdf(gammaDist3, beta));
    boost::math::gamma_distribution<> gammaDist4(prior_param[6], prior_param[7]);
    log_lhd += log(boost::math::pdf(gammaDist4, tilde_beta));
    log_lhd+= composite_log_likelihood_bc_tilde_woeps_infcentroids_inflated(clusters, centroids, t_end, a0, a1, beta, tilde_beta);
    return exp(log_lhd);
}