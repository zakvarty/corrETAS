#ifndef BAYESIANETAS_ADD_H
#define BAYESIANETAS_ADD_H

void rnormCpp(std::vector<double> &samples, int n, double mu, double sigma);

void discrete_sampleCpp(int size,int n, std::vector<double> &probs, std::vector<int> &samples);

int    rdiscrete(std::vector<double> &probs);

double dpoislog(double k, double lambda);

double var(std::vector<double> x);

double sd(std::vector<double> x);
#endif