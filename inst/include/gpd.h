#ifndef BAYESIANETAS_GPD_H
#define BAYESIANETAS_GPD_H 

double dgpdlog(double x, double nu, double xi, double mu);

double gpdllh(std::vector<double> &xs, double nu, double xi, double mu);

double gpdSurvivor(double x, double nu, double xi, double mu);

#endif