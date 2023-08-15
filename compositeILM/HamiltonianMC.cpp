//
// Created by yirao zhang on 2023-08-08.
//
#include "HamiltonianMC.h"
using namespace std;
vector<double> diff_log_target(Individual pop[], double a0, double a1, double beta){
    // domain space vector
    size_t n = 3;               // number of domain space variables
    vector<CppAD::AD<double>> ax(n); // vector of domain space variables
    ax[0] = a0;
    ax[1] = a1;
    ax[2] = beta;
    // declare independent variables and start recording operation sequence
    CppAD::Independent(ax);
    // range space vector
    size_t m = 1;               // number of ranges space variables
    vector< CppAD::AD<double> > ay(m); // vector of ranges space variables
    ay[0] = -log(post_dist_dens_logprior(pop, a0, a1, beta));     // record operations that compute ay[0]
    // store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<double> f(ax, ay);
    // compute derivative using operation sequence stored in f
    vector<double> jac(m * n); // Jacobian of f (m by n matrix)
    vector<double> x(n);       // domain space vector
    x[0] = a0;// argument value for computing derivative
    x[1] = a1;
    x[2] = beta;
    jac  = f.Jacobian(x);      // Jacobian for operation sequence
    return jac;
}
void leapfrog(Individual pop[], vector<double>& z, vector<double>& r, double e, int L) {//z-parameters(a0, a1, beta), r-anxillary variable, e-step size, L-number of iterations
    //double E = -log(post_dist_dens(pop, z[0], z[1], z[3])); //momentum
    vector<double> znxt(3, 0);
    vector<double> rnxt(3, 0);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<double> dist(0, 1);
    double rand = dist(gen);
    if (rand > 0.5) {
        e = -e;
    }
    vector<double> grad = diff_log_target(pop, z[0], z[1], z[2]);
    for (size_t i = 0; i < 3; ++i) {
        rnxt[i] = r[i] - 0.5 * e * grad[i];
    }
    znxt = z;
    for (size_t i = 0; i < L; i++) {
        //update znxt
        for (size_t i = 0; i < 3; ++i) {
            znxt[i] = znxt[i] + e * rnxt[i];
        }
        vector<double> grad = diff_log_target(pop, znxt[0], znxt[1], znxt[2]);
        //update rnxt
        for (size_t i = 0; i < 3; ++i) {
            rnxt[i] = rnxt[i] - e * grad[i];
        }
    }
    for (size_t i = 0; i < 3; ++i) {
        r[i] = rnxt[i] + 0.5 * e * grad[i];
        z[i] = znxt[i];
    }
}

void HMC(Individual pop[], double samples[][3]){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> dis(0,1);
    vector<double> r(3,0.1);
    vector<double> z(3, 0.1);
    for(size_t i=0; i<50000; i++){
        r[0] = dis(gen);
        r[1] = dis(gen);
        r[2] = dis(gen);
        vector<double> z0 = z;
        double H0 = -log(post_dist_dens_logprior(pop, z[0], z[1], z[2]))+0.5*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        leapfrog(pop, z, r, 0.1, 10);
        //cout <<z[0] <<" "<<z[1]<<" "<<z[2]<<endl;
        double Hn = -log(post_dist_dens_logprior(pop, z[0], z[1], z[2]))+0.5*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        double e = (double)(rand()%100001)/100000;
        if(e> exp(-Hn-(-H0))){
            z = z0;
        }
        samples[i][0] = exp(z[0]);
        samples[i][1] = exp(z[1]);
        samples[i][2] = exp(z[2]);
        cout << "e " <<e<< " "<<samples[i][0] <<" "<<samples[i][1]<<" "<<samples[i][2]<<endl;
    }
}