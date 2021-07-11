/*#include "randutils.hpp"*/
#include <R.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <Rmath.h>
//#include <add.h>
//#include <gpd.h>

using namespace std;

// -----------------------------------------------------------------------------
// ---------------- Function definitions: all     ------------------------------
// -----------------------------------------------------------------------------
// General stats functions: Defined in add.h and gpd.h
void rnormCpp(std::vector<double> &samples, int n, double mu, double sigma);

void discrete_sampleCpp(int size,int n, std::vector<double> &probs, std::vector<int> &samples);

int   rdiscrete(std::vector<double> &probs);

double dpoislog(double k, double lambda);

double var(std::vector<double> x);

double sd(std::vector<double> x);

// GPD related functions
double pgpd(double q, double nu, double xi, double mu);

double dgpdlog(double x, double nu, double xi, double mu);

double gpdllh(std::vector<double> &xs, double nu, double xi, double mu);

double gpdSurvivor(double x, double nu, double xi, double mu);

double gpd_to_gauss(double x_gpd, double nu, double xi, double mu);

double dnormlog(double x, double mu, double sd);

double dbvnlog(double x1, double x2, double mu1, double mu2, double sd1, double sd2, double rho);

double corrmag_weight(int p, int i, std::vector<double> &marks, std::vector<int> &branching, double mag_thresh, double numb, double ximb, double numa, double xima, double rho);

// Branching functions (cGPD)
void sampleBranchingcGPD(std::vector<double> &ts, std::vector<double> &marks, double mu, double logC, double alpha, double nu, double xi, std::vector<int> &branching);

void sampleBranchingcGPDdualmag(std::vector<double> &ts, std::vector<double> &marks,  double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double M0, std::vector<int> &branching);

void sampleBranchingcGPDcorrmag(std::vector<double> &ts, std::vector<double> &marks,  double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double rho, double M0, std::vector<int> &branching);

// ETAS functions

double hConditionalPosterior(std::vector<double> &ts, std::vector<double> &marks, std::vector<double> &zs, double maxTime, std::vector<double> &kappaevals, double nu, double xi, int constrainOmori);

double kappaConditonalPosterior(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &numtriggered, std::vector<double> &Hevals, double logC, double alpha);

double rhoConditionalPosterior(std::vector<double> &markstriggered, std::vector<double> &parent_marks, std::vector<double> &parent_distribution, double numb, double ximb, double numa, double xima, double rho, double M0);

bool is_psi_invalid(std::vector<double> marksbackground, std::vector<double> markstriggered, double numb, double ximb, double numa, double xima, double rho, double M0);

double psiConditionalPosterior(std::vector<double> marksbackground, std::vector<double> markstriggered, std::vector<double> parent_marks, std::vector<int> parent_distribution, double numb, double ximb, double numa, double xima, double rho, double M0);

void estimateETAS_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double num, double xim, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &nums, std::vector<double> &xims, std::vector<int> &Bs);

void estimateETASdual_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &numbs, std::vector<double> &ximbs,std::vector<double> &numas, std::vector<double> &ximas, std::vector<int> &Bs);

void estimateETAScorr_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double rho, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &numbs, std::vector<double> &ximbs,std::vector<double> &numas, std::vector<double> &ximas, std::vector<double> &rhos, std::vector<int> &Bs);

/* ---------------------------------------------------------------------------*/

// -----------------------------------------------------------------------------
// ---------------- General stats functions     --------------------------------
// -----------------------------------------------------------------------------
void rnormCpp(std::vector<double> &samples, int n, double mu, double sigma) {
	samples.clear();
	samples.reserve(n);

	GetRNGstate();
	/* double rnorm(double mu, double sigma); */
	for (int i = 0 ; i < n ; i ++) {
		double draw = rnorm(mu, sigma);
		samples.push_back(draw);
		Rprintf("\t samples[%d]:  %f \n", i, samples[i]);
	}
	PutRNGstate();
}


void discrete_sampleCpp(int size, int n, std::vector<double> &probs,  std::vector<int> &samples) {
	samples.clear();
	samples.reserve(size);

	/* Make vector of cumulative probabilities */
	std::vector<double> probsCum = probs;
	for(int i = 1 ; i < n ; i++) {
		probsCum[i] += probsCum[i-1] ;
	}

	GetRNGstate();
	for (int i = 0; i < size ; i++){
		/* take one sample */
		double U = runif(0,1);
		Rprintf("U[%d]:  %f \n", i, U);

		int sample = 0;
		for(int j = 1; j < n ; j++){
			if(probsCum[j] < U){
				sample += 1;
			}
		}
		/* append to samples */
		samples.push_back(sample);
		Rprintf("samples[%d]:  %f \n", i, samples[i]);
	}
	PutRNGstate();
}


int rdiscrete(std::vector<double> &probs) {
	GetRNGstate();
	int n = probs.size();
	double total = 0;
	for( int i = 0; i < n; i++){
		total += probs[i];
	}

	/* Make vector of cumulative probabilities */
	std::vector<double> probsCum = probs;
	probsCum[0] = 0;
	for(int i = 1 ; i < n ; i++) {
		probsCum[i] = probsCum[i-1] + probs[i-1] ;
	}
	/* Convert U to index of probs */
	double U = runif(0,total);
	int sample = 0;
	for(int j = 1; j < n ; j++){
		if(probsCum[j] < U){
			sample += 1;
		}
	}

	PutRNGstate();
	return(sample);
}


double dpoislog(double k, double lambda) {
	if (k==0) {
		return(k * log(lambda) - 0 - lambda);
	}
	if (k==1) {
		return(k * log(lambda) - 0 - lambda);
	}
	if (k==2) {
		return(k * log(lambda) - 0.6931472 - lambda);
	}
	if (k==3) {
		return(k * log(lambda) - 1.791759 - lambda);
	}
    return k * log(lambda) - lgamma(k + 1.0) - lambda;
}


double var(std::vector<double> x) {
     int size = x.size();

     double variance = 0;
     double mean = x[0];
     for (int i = 1; i < size; i++)
     {
          mean += x[i];
     }
	 mean = mean/size;
	 for (int i = 0; i < size; i++) {
		 variance += pow((x[i] - mean),2);
	 }

     return variance / (size - 1);
}


double sd(std::vector<double> x) {
	return sqrt(var(x));
}


// -----------------------------------------------------------------------------
// ---------------- GPD  related functions     ---------------------------------
// -----------------------------------------------------------------------------
double pgpd(double q, double nu, double xi, double mu) {
	double sig = nu/(1 + xi);
	double uep = mu - sig/xi;
	double p;

	if (q <= mu) {
		p = 0.0;
	} else if ( (xi < 0) && (q >= uep)) {
		p = 1;
	} else if ( xi == 0.0) {
		p = (1 - exp(-(q - mu) / nu));
	} else {
		p = (1 - pow(1 + (xi * (q - mu))/sig, -1/xi));
	}
	return(p);
}


double dgpdlog(double x, double nu, double xi, double mu) {
	if(x < mu){
		return(-9999999); //should be -Inf
	}
	if(xi==0) {
		return(-log(nu)-(x-mu)*(1+xi)/nu);
	}
	if(xi < 0) {
		if(x > (mu - nu/(xi*(1 + xi)))){
			return(-9999999);
		}
	}
	return(-log(nu) + log(1+xi) - (1/xi + 1) * log(1 + xi*(1+xi)*(x-mu)/nu));
}


double gpdllh(std::vector<double> &xs, double nu, double xi, double mu){
	int size = xs.size();
	if(nu / (1.0 + xi) < 0.0){
		return(-9999999);
	}
	double llh = 0;
	for (int i = 0; i < size ; i ++){
		llh += dgpdlog(xs[i],nu,xi,mu);
	}
	return(llh);
}


double gpdSurvivor(double x, double nu, double xi, double mu) {
	if(x<mu) {
		return(1);
	}
	if(xi==0) {
		return( exp(-(x-mu)*(1+xi)/nu) );
	}
	if(xi<0){
		if(x > (mu - nu/(xi * (1+xi)))) {
			return(0);
		}
	}
	return(pow((1 + (x - mu)*xi*(1+xi)/nu),(-1/xi)));
}


double gpd_to_gauss(double x_gpd, double nu, double xi, double mu){
	double x_gauss = 0.0;
    double x_unif;
    double x_gauss_approx;
	double epsilon;
	double sig = nu / (1 + xi);                // GPD scale parameter

	x_unif = pgpd(x_gpd, nu, xi, mu);
	if (x_unif < 0.999) {
		x_gauss = qnorm( /*p=*/ x_unif, /*mu=*/0.0, /*sigma=*/ 1.0, /*lower_tail=*/ true,/*log_p=*/false);
	} else {
		if (xi != 0) {
			x_gauss_approx = sqrt( (2.0/xi) * log( 1 + (xi/sig)*(x_gpd - mu)) );
			epsilon = (-log((4*M_PI/xi)*log(1+(xi/sig)*(x_gpd - mu)))) / ((4/xi)*log(1 + (xi/sig)*(x_gpd - mu)));
		} else {
			x_gauss_approx = sqrt(2*(x_gauss - mu)/sig);
			epsilon = log(2*M_PI*(x_gpd - mu)/sig) / (4*(x_gpd - mu)/sig);
		}
		x_gauss = x_gauss_approx + epsilon;
	}
	return(x_gauss);
}

double dnormlog(double x, double mu, double sd){
	double log_density = -0.5 * log(2 * M_PI) - log(sd) - 0.5 * pow(sd,-2) * pow(x - mu, 2); //final * was a +
	return(log_density);
}

double dbvnlog(double x1, double x2, double mu1, double mu2, double sd1, double sd2, double rho){

	double mu3 = mu2 + (sd2 / sd1) * rho * (x1 - mu1);
	double sd3 = sd2 * sqrt(1 - pow(rho, 2));

	double log_density = -log(2 * M_PI) ;
	log_density -= log(sd1); //this was missing but does not depend on rho
	log_density -= log(sd3); //this was missing and depends on rho
	log_density +=   (- 1.0 / (2 * pow(sd1, 2))) * pow(x1 - mu1, 2) ;
	log_density +=   (- 1.0 / (2 * pow(sd3, 2))) * pow(x2 - mu3, 2) ;

	return(log_density);
}

double corrmag_weight(int p, int i, std::vector<double> &marks, std::vector<int> &branching, double mag_thresh, double numb, double ximb, double numa, double xima, double rho){

	// NB: p is parent process number.  p \in 0:i-1
	//    respective parent mags is in marks[NA,0,1,...,i-2]
	//     i is the ordinal of current event. i in 1:n. (E.g. for first event i = 1, not 0.)
	//    respective mags are in marks[0,1,2,...,n-1]


	int n = marks.size();
	int hasParent;
	int parTrig;

	double nu_i;
	double xi_i;
	double nu_p;
	double xi_p;
	double nu_c = numa;
	double xi_c = xima;
	double log_density;
	double m_i;
	double z_i;
	double m_p;
	double z_p;


	std::vector<int> C;
	std::vector<double> m_c;
	std::vector<double> z_c;
	C.reserve(n);


	// Get magnitude parameters of current event
	hasParent =  (p > 0);
	nu_i = numa * hasParent + numb * (1 - hasParent);
	xi_i = xima * hasParent + ximb * (1 - hasParent);

	// Get magnitude parameters of parent event
	if(hasParent){
	    parTrig = (branching[p-1] > 0);
	    nu_p = numa * parTrig + numb * (1 - parTrig);
	    xi_p = xima * parTrig + ximb * (1 - parTrig);
	}

	// Get locations of child events in branching
	for(int j = i; j < n; j++){
		if(branching[j] == i){
			C.push_back(j);
		}
	}

	int modC = C.size();
	int hasChildren = (modC > 0);

	// Get child magnitude(s) on GPD and Gaussian margins
	if(hasChildren){
		m_c.reserve(modC);
		z_c.reserve(modC);

		m_c.clear();
		z_c.clear();

		for(int j = 0; j < modC; j++){
			m_c.push_back(marks[C[j]]);
			z_c.push_back(gpd_to_gauss(m_c[j], numa, xima, mag_thresh)); //???
		}
	}
	// Get current magnitude on GPD and Gaussian margins
	m_i = marks[i-1];
	z_i = gpd_to_gauss(m_i, nu_i, xi_i, mag_thresh);
	// Get parent magnitude on GPD and Gaussian margins
	if(hasParent){
		m_p = marks[p-1];
		z_p = gpd_to_gauss(m_p, nu_p, xi_p, mag_thresh);
	}

	// Calculate log-density
	/*log_density = dgpdlog(m_i, nu_i, xi_i, mag_thresh);

	if(hasParent){
		log_density += dgpdlog(m_p, nu_p, xi_p, mag_thresh);
		log_density += dbvnlog(z_p, z_i, 0.0, 0.0, 1.0, 1.0, rho);

	}

	if(hasChildren){
		for(int j = 0; j < modC; j++){
			log_density += dgpdlog(m_c[j], nu_c, xi_c, mag_thresh);
			log_density += dbvnlog(z_i, z_c[j], 0.0, 0.0, 1.0, 1.0, rho);
		}
	} */

	log_density = 0.0;

	if(hasParent){
		log_density += dnormlog(z_i, rho * z_p, sqrt(1 - pow(rho,2))) ; 
		log_density -= dnormlog(z_i, 0 ,1);
	} 

	log_density += dgpdlog(m_i, nu_i, xi_i, mag_thresh);

	if(hasChildren){
		for(int j = 0; j < modC; j++){
			log_density += dnormlog(z_c[j], rho * z_i, sqrt(1 - pow(rho,2)));
			log_density -= dnormlog(z_c[j], 0, 1);
			log_density += dgpdlog(m_c[j], nu_c, xi_c, mag_thresh);
		}
	}

	return(exp(log_density));
}


// -----------------------------------------------------------------------------
// ---------------- Branching functions        ---------------------------------
// -----------------------------------------------------------------------------


/******** BRANCHING FUNCTIONS cGPD **/
void sampleBranchingcGPD(std::vector<double> &ts, std::vector<double> &marks,  double mu, double logC, double alpha, double nu, double xi, std::vector<int> &branching) {
	int n = ts.size();
	int parent;
	double temp;
	double C = exp(logC);
	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	GetRNGstate();

	branching.clear();
	branching.reserve(n);
	branching.push_back(0);

	std::vector<double> lambdas;
	lambdas.reserve(n);

	std::vector<double> kappaevals;
	kappaevals.reserve(n);

	//precompute
	for (int i = 0; i < n; i++) {
		kappaevals.push_back(C*exp(alpha*(marks[i]-Mbar)));
	}

	for (int i = 1; i < n; i++) {
		lambdas.clear();
		lambdas.push_back(mu);

		for (int j = 0; j < i; j++) { //check this iteratures through...
			temp = kappaevals[j]*exp(dgpdlog(ts[i]-ts[j], nu, xi, 0));
			lambdas.push_back(temp);
		}

		parent = rdiscrete(lambdas);
		branching.push_back(parent);
	}
	PutRNGstate();
}

void sampleBranchingcGPDdualmag(std::vector<double> &ts, std::vector<double> &marks,  double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double M0, std::vector<int> &branching) {
	int n = ts.size();
	int parent;
	double temp;
	double C = exp(logC);
	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	GetRNGstate();

	branching.clear();
	branching.reserve(n);
	branching.push_back(0);

	std::vector<double> weights;
	weights.reserve(n);

	std::vector<double> lambdas;
	lambdas.reserve(n);

	std::vector<double> weightedlambdas;
	weightedlambdas.reserve(n);

	std::vector<double> kappaevals;
	kappaevals.reserve(n);

	//precompute
	for (int i = 0; i < n; i++) {
		kappaevals.push_back(C*exp(alpha*(marks[i]-Mbar)));
	}

	for (int i = 1; i < n; i++) {                                                   /* for the second, third,... last event (i)*/
		weights.clear();
		lambdas.clear();
		weightedlambdas.clear();

		for (int j = 0; j < (i+1); j ++) {												/*for each possible parent process  j = 0 ,... i-1 */
			int hasParent = std::min(1, j);										/*Is event i a triggered event with this parent process?*/
			double num = numa*hasParent + numb*(1-hasParent);					/* get the triggered/background parameters as required */
			double xim = xima*hasParent + numb*(1-hasParent);

			temp = exp(dgpdlog(marks[i], num, xim, M0));					        /* convert this to a weight by looking up the density of magnitude 2,...,n and multiplying by 2,3,... n */
			weights.push_back(temp);	 											/* record this weight */
		}

		lambdas.push_back(mu);														/* for parent process 0  */
		for (int j = 0; j < i; j++) {                                           	/* and for parent processes j = 1,...,i-1 */
			temp = kappaevals[j]*exp(dgpdlog(ts[i]-ts[j], nut, xit, 0));        /* calculate \lambda_j at t_{i+1} (stored in ts[i]) */
			lambdas.push_back(temp);                                                /* NB:  \kappa(m_j) is stored in kappaevals[j-1] */
		}

		for (int j = 0; j < (i+1); j++) {                                           /* zero indexing is a nightmare. */
			weightedlambdas.push_back(weights[j]*lambdas[j]);
		}

		parent = rdiscrete(weightedlambdas);
		branching.push_back(parent);
	}
	PutRNGstate();
}


void sampleBranchingcGPDcorrmag(std::vector<double> &ts, std::vector<double> &marks,  double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double rho, double M0, std::vector<int> &branching) {
	int n = ts.size();
	int parent;
	double temp;
	double C = exp(logC);
	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	GetRNGstate();

	branching.clear();
	branching.reserve(n);
	branching.push_back(0);

	std::vector<double> weights;
	weights.reserve(n);

	std::vector<double> lambdas;
	lambdas.reserve(n);

	std::vector<double> weighted_lambdas;
	weighted_lambdas.reserve(n);

	std::vector<double> kappaevals;
	kappaevals.reserve(n);


	//precompute
	for (int i = 0; i < n; i++) {
		kappaevals.push_back(C*exp(alpha*(marks[i]-Mbar)));
	}


	for (int i = 1; i < n; i++) {                                                   /* for the second, third,... last event (i)*/
		weights.clear();
		lambdas.clear();
		weighted_lambdas.clear();

		for (int j = 0; j < (i + 1); j ++){ // ??
         //corrmag_weight(p (parent process #), i (1-index event #),  marks, branching, mag_thresh,  numb, double ximb, double numa, double xima, double rho){
			temp = corrmag_weight(j, i+1, marks, branching, M0, numb, ximb, numa, xima, rho);
			weights.push_back(temp);
		}


		lambdas.push_back(mu);														/* for parent process 0  */
		for (int j = 0; j < i; j++) {                                           	/* and for parent processes j = 1,...,i-1 */
			temp = kappaevals[j]*exp(dgpdlog(ts[i]-ts[j], nut, xit, 0));        /* calculate \lambda_j at t_{i+1} (stored in ts[i]) */
			lambdas.push_back(temp);                                                /* NB:  \kappa(m_j) is stored in kappaevals[j-1] */
		}

		for (int j = 0; j < (i+1); j++) {                                           /* zero indexing is a nightmare. */
			weighted_lambdas.push_back(weights[j]*lambdas[j]);
		}

		parent = rdiscrete(weighted_lambdas);
		branching.push_back(parent);
	}
	PutRNGstate();
}

// -----------------------------------------------------------------------------
// ---------------- ETAS functions        ---------------------------------
// -----------------------------------------------------------------------------

double hConditonalPosterior(std::vector<double> &ts, std::vector<double> &marks, std::vector<double> &z, double maxTime, std::vector<double> &kappaevals, double nu, double xi, int constrainOmori) {
	if (nu/(1+xi) <= 0) {
		return(-9999999);
	}
	if (constrainOmori == 1) {
		if (xi <= 0 || nu <= 0){
			return(-9999999);
		}
	}

	double loglik = 0, temp;
	int n = ts.size();
	for (int i = 0; i < n; i++) {
		temp = 1 - gpdSurvivor(maxTime-ts[i], nu, xi, 0);
		loglik = loglik - kappaevals[i]*temp;
	}
	for (unsigned int i = 0; i < z.size(); i++) {
		loglik += dgpdlog(z[i], nu, xi, 0);
	}
	return(loglik);
}

double kappaConditionalPosterior(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &numtriggered, std::vector<double> &Hevals, double logC,  double alpha) {
	if (alpha < 0 || alpha > 10) {
		return(-9999999);
	}

	double loglik = 0;
	int n = marks.size();
	double temp;
	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	for (int i = 0 ; i < n; i++) {
		temp = exp(logC + alpha*(marks[i]-Mbar));
		temp = temp * Hevals[i];
		loglik += dpoislog(numtriggered[i], temp);
	}
	return(loglik);
}

double rhoConditionalPosterior(std::vector<double> markstriggered, std::vector<double> parent_marks, std::vector<int> parent_distribution, double numb, double ximb, double numa, double xima, double rho, double M0){
	if (rho < -1 || rho > 1){
		return(-9999999);
	}

	double loglik = 0;
	int n = markstriggered.size();

	for (int i = 0; i < n ; i++){
		double z_i = gpd_to_gauss(markstriggered[i], numa, xima, M0);

		int parTrig = parent_distribution[i] > 0;
		double nu_p =  parTrig * numa + (1 - parTrig) * numb;
		double xi_p =  parTrig * xima + (1 - parTrig) * ximb;
		double z_p = gpd_to_gauss(parent_marks[i], nu_p, xi_p, M0);

		//loglik += dnormlog(z_p, 0, 1);
		//loglik += dnormlog(z_i, rho * z_p, sqrt(1 - pow(rho,2)));
		loglik += dbvnlog(z_p, z_i, 0.0, 0.0, 1.0, 1.0, rho);
	}
	return(loglik);
}

// EDITING ZONE
bool is_psi_invalid(std::vector<double> marksbackground, std::vector<double> markstriggered, double numb, double ximb, double numa, double xima, double rho, double M0){
	// check valid correlation
	if (rho < -1 || rho > 1){
		return(TRUE);
	}

	// check valid background magnitude parameters
	if(numb / (1.0 + ximb) < 0.0){
		return(TRUE);
	}

	if(ximb < 0){
		double uppermb = (M0 - numb/(ximb*(1 + ximb)));
		for(int i = 0; i < marksbackground.size(); i ++){
			if(marksbackground[i] >= uppermb){
				return(TRUE);
			}
		}
	}

	// check valid aftershock magnitude parameters
	if(numa / (1.0 + xima) < 0.0){
		return(TRUE);
	}

	if(xima < 0){
		double upperma = (M0 - numa/(xima*(1 + xima)));
		for(int i = 0; i < markstriggered.size(); i ++){
			if(markstriggered[i] >= upperma){
				return(TRUE);
			}
		}
	}

	return(FALSE);
}

double psiConditionalPosterior(std::vector<double> marksbackground, std::vector<double> markstriggered, std::vector<double> parent_marks, std::vector<int> parent_distribution, double numb, double ximb, double numa, double xima, double rho, double M0){
	
	if(is_psi_invalid(marksbackground, markstriggered, numb, ximb, numa, xima, rho, M0)){
		return(-99999999);
	}
	
	double loglik = 0;

	loglik += gpdllh(marksbackground, numb, ximb, M0);
	loglik += gpdllh(markstriggered, numa, xima, M0);

	for (int i = 0; i < markstriggered.size() ; i++){
		double z_i = gpd_to_gauss(markstriggered[i], numa, xima, M0);

		int parTrig = parent_distribution[i] > 0;
		double nu_p =  parTrig * numa + (1 - parTrig) * numb;
		double xi_p =  parTrig * xima + (1 - parTrig) * ximb;
		double z_p = gpd_to_gauss(parent_marks[i], nu_p, xi_p, M0);

		double term1 = dnormlog(z_i, rho * z_p, sqrt(1.0 - pow(rho,2)));
		double term2 = dnormlog(z_i, 0.0, 1.0);
		
		loglik += term1; 
		loglik -= term2;
	}
	return(loglik);
}
// EDITING ZONE

void estimateETAS_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double num, double xim, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &nums, std::vector<double> &xims, std::vector<int> &Bs){

	int n = ts.size(), numbackground;
	double currposterior,newposterior, newnut,newxit,newlogC,newalpha,newnum,newxim;

	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	double mualpha = muprior[0];
	double mubeta  = muprior[1]; //prior params

	std::vector<double> kappaevals;
	kappaevals.reserve(n);
	std::vector<double> Hevals;
	Hevals.reserve(n);

    GetRNGstate();

	Rprintf("Using the following prior for mu: \n");
	Rprintf("\t mu ~ Gamma(alpha = %f, beta = %f) \n", mualpha, mubeta);

	Rprintf("Metropolis steps will have the following standard deviations: \n");
    Rprintf("\t logC:  %f \n", etasstepSds[0]);
	Rprintf("\t alpha: %f \n", etasstepSds[1]);
	Rprintf("\t nu_t:     %f \n", etasstepSds[2]);
	Rprintf("\t xi_t:     %f \n", etasstepSds[3]);
	Rprintf("\n");
	Rprintf("\t nu_m:     %f \n", magstepSds[0]);
	Rprintf("\t xi_m:     %f \n", magstepSds[1]);



	for (int s=0; s<sims; s++) {


		if (Bfixed != 1) {
		sampleBranchingcGPD(ts,marks,mu,logC,alpha,nut,xit,branching);
		}
		if (Bsamples == 1) {
			for (int i=0; i<n; i++) {
			Bs.push_back(branching[i]);
			}
		}


		std::vector<int> numtriggered(n,0);
		std::vector<double> z;
		z.reserve(n);

		numbackground = 0;
		for (int i=0; i<n; i++) {
			if (branching[i] > 0) {
				numtriggered[branching[i]-1]++;
				z.push_back(ts[i]-ts[branching[i]-1]);
			} else {
				numbackground++;
			}
		}

		mu = rgamma(mualpha+numbackground, 1/(mubeta+maxTime));
		mus.push_back(mu);

		kappaevals.clear();
		for (int i=0; i<n; i++) {
			kappaevals.push_back(exp(logC + alpha*(marks[i]-Mbar)));
		}

		currposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,nut,xit,constrainOmori);
		for (int i=0; i < nuxiSteps; i++) {
			newnut = nut + etasstepSds[2]*rnorm(0,1);
			newxit = xit + etasstepSds[3]*rnorm(0,1);
			newposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,newnut,newxit,constrainOmori);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				nut = newnut;
				xit = newxit;
				currposterior = newposterior;
			}
		}
		nuts.push_back(nut);
		xits.push_back(xit);

		Hevals.clear();
		for (int i=0; i<n; i++) {
			Hevals.push_back(1 - gpdSurvivor(maxTime-ts[i], nut, xit, 0));
		}

		currposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,logC,alpha);
		for (int i=0; i < CaSteps; i++) {
			newlogC = logC + etasstepSds[0]*rnorm(0,1);
			newalpha = alpha + etasstepSds[1]*rnorm(0,1);
			newposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,newlogC,newalpha);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				logC = newlogC;
				alpha = newalpha;
				currposterior = newposterior;
			}
		}
		logCs.push_back(logC);
		alphas.push_back(alpha);

		currposterior = gpdllh(marks,num,xim,M0);
		for (int i=0; i < magSteps; i++) {
			newnum = num + magstepSds[0]*rnorm(0,1);
			newxim = xim + magstepSds[1]*rnorm(0,1);
			newposterior = gpdllh(marks,newnum,newxim,M0);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				num = newnum;
				xim = newxim;
				currposterior = newposterior;
			}
		}
		nums.push_back(num);
		xims.push_back(xim);

		if (s % 100 == 0) {
			Rprintf("Generated %d samples...\n",s);
		}
	}
	PutRNGstate();
	Rprintf("Done.");
}


void estimateETASdual_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &numbs, std::vector<double> &ximbs,std::vector<double> &numas, std::vector<double> &ximas, std::vector<int> &Bs){

	int n = ts.size(), numbackground;
	double currposterior,newposterior, newnut,newxit,newlogC,newalpha,newnumb,newximb,newnuma,newxima;

	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	double mualpha = muprior[0];
	double mubeta  = muprior[1]; //prior params of gamma on mu

	std::vector<double> kappaevals;
	kappaevals.reserve(n);
	std::vector<double> Hevals;
	Hevals.reserve(n);

    GetRNGstate();

	Rprintf("Using the following prior for mu: \n");
	Rprintf("\t mu ~ Gamma(alpha = %f, beta = %f) \n", mualpha, mubeta);

	Rprintf("Metropolis steps will have the following standard deviations: \n");
    Rprintf("\t logC:  %f \n", etasstepSds[0]);
	Rprintf("\t alpha: %f \n", etasstepSds[1]);
	Rprintf("\t nu_t:     %f \n", etasstepSds[2]);
	Rprintf("\t xi_t:     %f \n", etasstepSds[3]);
	Rprintf("\n");
	Rprintf("\t nu_mb:     %f \n", magstepSds[0]);
	Rprintf("\t xi_mb:     %f \n", magstepSds[1]);
	Rprintf("\t nu_ma:     %f \n", magstepSds[2]);
	Rprintf("\t xi_ma:     %f \n", magstepSds[3]);


	for (int s=0; s<sims; s++) {

		if (Bfixed != 1) {
		sampleBranchingcGPDdualmag(ts,marks,mu,logC,alpha,nut,xit,numb,ximb,numa,xima,M0,branching);
		}
		if (Bsamples == 1) {
			for (int i=0; i<n; i++) {
			Bs.push_back(branching[i]);
			}
		}

		std::vector<int> numtriggered(n,0);
		std::vector<double> z;
		std::vector<double> marksbackground;
		std::vector<double> markstriggered;
		z.reserve(n);
		marksbackground.reserve(n);
		markstriggered.reserve(n);

		numbackground = 0;
		for (int i=0; i<n; i++) {
			if (branching[i] > 0) {
				numtriggered[branching[i]-1]++;
				z.push_back(ts[i]-ts[branching[i]-1]);
				markstriggered.push_back(marks[i]);
			} else {
				numbackground++;
				marksbackground.push_back(marks[i]);
			}
		}

		mu = rgamma(mualpha+numbackground, 1/(mubeta+maxTime));
		mus.push_back(mu);

		kappaevals.clear();
		for (int i=0; i<n; i++) {
			kappaevals.push_back(exp(logC + alpha*(marks[i]-Mbar)));
		}

		currposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,nut,xit,constrainOmori);
		for (int i=0; i < nuxiSteps; i++) {
			newnut = nut + etasstepSds[2]*rnorm(0,1);
			newxit = xit + etasstepSds[3]*rnorm(0,1);
			newposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,newnut,newxit,constrainOmori);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				nut = newnut;
				xit = newxit;
				currposterior = newposterior;
			}
		}
		nuts.push_back(nut);
		xits.push_back(xit);

		Hevals.clear();
		for (int i=0; i<n; i++) {
			Hevals.push_back(1 - gpdSurvivor(maxTime-ts[i], nut, xit, 0));
		}

		currposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,logC,alpha);
		for (int i=0; i < CaSteps; i++) {
			newlogC = logC + etasstepSds[0]*rnorm(0,1);
			newalpha = alpha + etasstepSds[1]*rnorm(0,1);
			newposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,newlogC,newalpha);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				logC = newlogC;
				alpha = newalpha;
				currposterior = newposterior;
			}
		}
		logCs.push_back(logC);
		alphas.push_back(alpha);

		Rprintf("\t s:     %d ", s);
    Rprintf("\t n_back:     %d ", numbackground);
    Rprintf("\t n_trig:     %d  \n", n - numbackground);

		currposterior = gpdllh(marksbackground,numb,ximb,M0);
		for (int i=0; i < magSteps; i++) {
			newnumb = numb + magstepSds[0]*rnorm(0,1);
			newximb = ximb + magstepSds[1]*rnorm(0,1);
			newposterior = gpdllh(marksbackground,newnumb,newximb,M0);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				numb = newnumb;
				ximb = newximb;
				currposterior = newposterior;
			}
		}
		numbs.push_back(numb);
		ximbs.push_back(ximb);

		currposterior = gpdllh(markstriggered,numa,xima,M0);
		for (int i=0; i < magSteps; i++) {
			newnuma = numa + magstepSds[2]*rnorm(0,1);
			newxima = xima + magstepSds[3]*rnorm(0,1);
			newposterior = gpdllh(markstriggered,newnuma,newxima,M0);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				numa = newnuma;
				xima = newxima;
				currposterior = newposterior;
			}
		}
		numas.push_back(numa);
		ximas.push_back(xima);


		if (s % 100 == 0) {
			Rprintf("Generated %d samples...\n",s);
		}
	}
	PutRNGstate();
	Rprintf("Done.");
}


void estimateETAScorr_Cpp(std::vector<double> &ts, std::vector<double> &marks, std::vector<int> &branching, double M0, double maxTime, int sims, int Bsamples, int Bfixed, int nuxiSteps, int CaSteps, int magSteps, double mu, double logC, double alpha, double nut, double xit, double numb, double ximb, double numa, double xima, double rho, std::vector<double> &etasstepSds, std::vector<double> &magstepSds, std::vector<double> &muprior, std::vector<double> &logCprior, int constrainOmori, std::vector<double> &mus, std::vector<double> &logCs, std::vector<double> &alphas, std::vector<double> &nuts, std::vector<double> &xits, std::vector<double> &numbs, std::vector<double> &ximbs,std::vector<double> &numas, std::vector<double> &ximas, std::vector<double> &rhos, std::vector<int> &Bs){

	int n = ts.size(), numbackground, check_count;
	double currposterior,newposterior,newnut,newxit,newlogC,newalpha,newnumb,newximb,newnuma,newxima,newrho;

	double Mbar = accumulate(marks.begin(),marks.end(),0.0)/((double)n);

	double mualpha = muprior[0];
	double mubeta  = muprior[1]; //prior params of gamma on mu
	double logC_mean = logCprior[0];
	double logC_sd   = logCprior[1]; //prior params of normal on logC

	std::vector<double> kappaevals;
	kappaevals.reserve(n);
	std::vector<double> Hevals;
	Hevals.reserve(n);

    GetRNGstate();

	Rprintf("Using the following prior for mu: \n");
	Rprintf("\t mu ~ Gamma(alpha = %f, beta = %f) \n", mualpha, mubeta);

	Rprintf("Initial C++ value for rho: %f \n", rho);

	Rprintf("Metropolis steps will have the following standard deviations: \n");
    Rprintf("\t logC:  %f \n", etasstepSds[0]);
	Rprintf("\t alpha: %f \n", etasstepSds[1]);
	Rprintf("\t nu_t:     %f \n", etasstepSds[2]);
	Rprintf("\t xi_t:     %f \n", etasstepSds[3]);
	Rprintf("\n");
	Rprintf("\t nu_mb:     %f \n", magstepSds[0]);
	Rprintf("\t xi_mb:     %f \n", magstepSds[1]);
	Rprintf("\t nu_ma:     %f \n", magstepSds[2]);
	Rprintf("\t xi_ma:     %f \n", magstepSds[3]);
    Rprintf("\t   rho:     %f \n", magstepSds[4]);

	for (int s=0; s<sims; s++) {
		//Rprintf("At iteration %d, location 1, rho = %f \n", s, rho);
		if (Bfixed != 1) {
		sampleBranchingcGPDcorrmag(ts,marks,mu,logC,alpha,nut,xit,numb,ximb,numa,xima,rho,M0,branching);
		}
		if (Bsamples == 1) {
			// resample vector
			for (int i=0; i<n; i++) {
			Bs.push_back(branching[i]);
			}
			// Check not all zero
			check_count = 0;
			for (int i=0; i<n; i++) {
			check_count += branching[i];
			}

			int resample_count = 0;
			// Repeat until at least one triggered
			while(check_count == 0){
				resample_count ++;
				Rprintf("B resampled %d times \n", resample_count);
				// resample vector
				for (int i=0; i<n; i++) {
					Bs.push_back(branching[i]);
				}
				// Check not all zero
				check_count = 0;
				for (int i=0; i<n; i++) {
					check_count += branching[i];
				}
			}
		}

		std::vector<int> numtriggered(n,0);
		std::vector<double> z;
		std::vector<double> marksbackground;
		std::vector<double> markstriggered;
		std::vector<double> parent_marks;
		std::vector<int> parent_distribution;
		z.reserve(n);
		marksbackground.reserve(n);
		markstriggered.reserve(n);
		parent_marks.reserve(n);
		parent_distribution.reserve(n);

		numbackground = 0;
		for (int i=0; i<n; i++) {
			if (branching[i] > 0) {
				numtriggered[branching[i]-1]++;
				z.push_back(ts[i]-ts[branching[i]-1]);
				markstriggered.push_back(marks[i]);
				parent_marks.push_back(marks[branching[i]-1]);
				parent_distribution.push_back(branching[branching[i]-1] > 0);
			} else {
				numbackground++;
				marksbackground.push_back(marks[i]);
			}
		}

		//Rprintf("not updating mu");
		mu = rgamma(mualpha+numbackground, 1/(mubeta+maxTime));
		mus.push_back(mu);

		kappaevals.clear();
		for (int i=0; i<n; i++) {
			kappaevals.push_back(exp(logC + alpha*(marks[i]-Mbar)));
		}

		currposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,nut,xit,constrainOmori);
		for (int i=0; i < nuxiSteps; i++) {
			newnut = nut + etasstepSds[2]*rnorm(0,1);
			newxit = xit + etasstepSds[3]*rnorm(0,1);
			newposterior = hConditonalPosterior(ts,marks,z,maxTime,kappaevals,newnut,newxit,constrainOmori);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				nut = newnut;
				xit = newxit;
				currposterior = newposterior;
			}
		}
		nuts.push_back(nut);
		xits.push_back(xit);

		Hevals.clear();
		for (int i=0; i<n; i++) {
			Hevals.push_back(1 - gpdSurvivor(maxTime-ts[i], nut, xit, 0));
		}

		currposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,logC,alpha);
		currposterior += dnormlog(logC, logC_mean, logC_sd);
		for (int i=0; i < CaSteps; i++) {
			newlogC = logC + etasstepSds[0]*rnorm(0,1);
			newalpha = alpha + etasstepSds[1]*rnorm(0,1);
			newposterior = kappaConditionalPosterior(ts,marks,numtriggered,Hevals,newlogC,newalpha);
			newposterior += dnormlog(newlogC, logC_mean, logC_sd);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				logC = newlogC;
				alpha = newalpha;
				currposterior = newposterior;
			}
		}
		logCs.push_back(logC);
		alphas.push_back(alpha);

		//currposterior = gpdllh(marksbackground,numb,ximb,M0);
		//for (int i=0; i < magSteps; i++) {
		//	newnumb = numb + magstepSds[0]*rnorm(0,1);
		//	newximb = ximb + magstepSds[1]*rnorm(0,1);
		//	newposterior = gpdllh(marksbackground,newnumb,newximb,M0);
		//	if (runif(0,1) < exp(newposterior-currposterior)) {
		//		numb = newnumb;
		//		ximb = newximb;
		//		currposterior = newposterior;
		//	}
		//}
		//numbs.push_back(numb);
		//ximbs.push_back(ximb);
		currposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, rho, M0);
		for (int i=0; i < magSteps; i++) {
			newnumb = numb + magstepSds[0]*rnorm(0,1);
			newximb = ximb + magstepSds[1]*rnorm(0,1);
			newposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, newnumb, newximb, numa, xima, rho, M0);
			if (runif(0,1) < exp(newposterior-currposterior)) {
				numb = newnumb;
				ximb = newximb;
				currposterior = newposterior;
			}
		}
		numbs.push_back(numb);
		ximbs.push_back(ximb);


		if(numbackground < n){
			//currposterior = gpdllh(markstriggered,numa,xima,M0);
			//for (int i=0; i < magSteps; i++) {
			//	newnuma = numa + magstepSds[2]*rnorm(0,1);
			//	newxima = xima + magstepSds[3]*rnorm(0,1);
			//	newposterior = gpdllh(markstriggered,newnuma,newxima,M0);
			//	if (runif(0,1) < exp(newposterior-currposterior)) {
			//		numa = newnuma;
			//		xima = newxima;
			//		currposterior = newposterior;
			//	}
			//}

			currposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, rho, M0);
			for (int i=0; i < magSteps; i++) {
				newnuma = numa + magstepSds[2]*rnorm(0,1);
				newxima = xima + magstepSds[3]*rnorm(0,1);
				newposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, numb, ximb, newnuma, newxima, rho, M0);
				if (runif(0,1) < exp(newposterior-currposterior)) {
					numa = newnuma;
					xima = newxima;
					currposterior = newposterior;
				}
			}
		}
		numas.push_back(numa);
		ximas.push_back(xima);

		//Rprintf("At iteration %d, location 5, rho = %f \n", s, rho);
		if(numbackground < n){
			//currposterior = rhoConditionalPosterior(markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, rho, M0);
			currposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, rho, M0);
			
			//Rprintf("At iteration %d, location 5.1, rho = %f \n", s, rho);
			for (int i=0; i < magSteps; i++) {
				double temp = magstepSds[4]*rnorm(0,1);
				newrho = rho + temp;
				//Rprintf("Current rho: %f \n step: %f \n Proposed rho: %f \n", rho, temp, newrho);
				//newposterior = rhoConditionalPosterior(markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, newrho, M0);
				newposterior = psiConditionalPosterior(marksbackground, markstriggered, parent_marks, parent_distribution, numb, ximb, numa, xima, newrho, M0);
				if (runif(0,1) < exp(newposterior-currposterior)) {
					rho = newrho;
					currposterior = newposterior;
					Rprintf("Accepted rho = %f \n", rho);
				}
			}
		}
		//Rprintf("At iteration %d, location 6: rho = %f \t n_back = %d, n_trig = %d\n", s, rho, numbackground, n - numbackground);
		Rprintf("At iteration %d: rho = %f, \t n_back = %d, \t n_trig = %d. \n", s, rho, numbackground, n - numbackground);
		rhos.push_back(rho);

		if (s % 100 == 0) {
			Rprintf("Generated %d samples...\n",s);
		}
	}

	//for (int i = 0 ; i < 20; i++){
	//	Rprintf("rhos[%s] = %f", i , rhos[i]);
	//}

	PutRNGstate();
	Rprintf("Done.");
}

extern "C" {
  // [[ export() ]]
  void estimateETAS_C(double *ts, double *marks, int *branching, int *n, double *M0, double *maxTime, int *sims, int *Bsamples, int *Bfixed, int *Bsize, int *nuxiSteps, int *CaSteps, int *magSteps, double *mu, double *logC, double *alpha, double *nut, double *xit, double *num, double *xim, double *etasstepSds, double *magstepSds, double *muprior, int *constrainOmori, double *mus, double *logCs, double *alphas, double *nuts, double *xits, double *nums, double *xims, int *Bs) {
		//Rprintf("Location (1) in estimateETAS_C");
		std::vector<double> vts(ts,ts + *n);
		std::vector<double> vmarks(marks,marks + *n);
		std::vector<int> vbranching(branching,branching + *n);
		std::vector<double> vetasstepSds(etasstepSds, etasstepSds + 4);
		std::vector<double> vmagstepSds(magstepSds, magstepSds + 2);
		std::vector<double> vmuprior(muprior, muprior + 2);
		//Rprintf("Location (2) in estimateETAS_C");
   	 	std::vector<double> vmus;
    	std::vector<double> vlogCs;
    	std::vector<double> valphas;
   	 	std::vector<double> vnuts;
    	std::vector<double> vxits;
		std::vector<double> vnums;
		std::vector<double> vxims;
		std::vector<int> vBs;
    	vmus.reserve(*sims);
    	vlogCs.reserve(*sims);
    	valphas.reserve(*sims);
    	vnuts.reserve(*sims);
    	vxits.reserve(*sims);
		vnums.reserve(*sims);
    	vxims.reserve(*sims);
		vBs.reserve(*Bsize);
		//Rprintf("Location (3) in estimateETAS_C");
		estimateETAS_Cpp(vts, vmarks, vbranching, *M0, *maxTime, *sims, *Bsamples, *Bfixed, *nuxiSteps, *CaSteps, *magSteps, *mu, *logC, *alpha, *nut, *xit, *num, *xim, vetasstepSds, vmagstepSds, vmuprior, *constrainOmori, vmus, vlogCs, valphas, vnuts, vxits, vnums, vxims, vBs);
		//Rprintf("Location (4) in estimateETAS_C");
		//Rprintf("Still need to check estimateETAS(C++)");
		//Rprintf("Location (5) in estimateETAS_C");
		std::copy(vmus.begin(),vmus.end(),mus);
        std::copy(vlogCs.begin(),vlogCs.end(),logCs);
        std::copy(valphas.begin(),valphas.end(),alphas);
		std::copy(vnuts.begin(),vnuts.end(),nuts);
        std::copy(vxits.begin(),vxits.end(),xits);
		std::copy(vnums.begin(),vnums.end(),nums);
        std::copy(vxims.begin(),vxims.end(),xims);
		std::copy(vBs.begin(), vBs.end(), Bs);
		//Rprintf("Location (6) in estimateETAS_C");
		return;
	}

  // [[ export() ]]
  void estimateETASdual_C(double *ts, double *marks, int *branching, int *n, double *M0, double *maxTime, int *sims, int *Bsamples, int *Bsize, int *Bfixed, int *nuxiSteps, int *CaSteps, int *magSteps,  double *mu, double *logC, double *alpha, double *nut, double *xit, double *numb, double *ximb, double *numa, double *xima, double *etasstepSds, double *magstepSds, double *muprior, int *constrainOmori, double *mus, double *logCs, double *alphas, double *nuts, double *xits, double *numbs, double *ximbs, double *numas, double *ximas,int *Bs){
		std::vector<double> vts(ts,ts + *n);
		std::vector<double> vmarks(marks,marks + *n);
		std::vector<int> vbranching(branching,branching + *n);
		std::vector<double> vetasstepSds(etasstepSds, etasstepSds + 4);
		std::vector<double> vmagstepSds(magstepSds, magstepSds + 4);
		std::vector<double> vmuprior(muprior, muprior + 2);

   	 	std::vector<double> vmus;
    	std::vector<double> vlogCs;
    	std::vector<double> valphas;
   	 	std::vector<double> vnuts;
    	std::vector<double> vxits;
		std::vector<double> vnumbs;
		std::vector<double> vximbs;
		std::vector<double> vnumas;
		std::vector<double> vximas;
		std::vector<int> vBs;
    	vmus.reserve(*sims);
    	vlogCs.reserve(*sims);
    	valphas.reserve(*sims);
    	vnuts.reserve(*sims);
    	vxits.reserve(*sims);
		vnumbs.reserve(*sims);
    	vximbs.reserve(*sims);
		vnumas.reserve(*sims);
    	vximas.reserve(*sims);
		vBs.reserve(*Bsize);


		estimateETASdual_Cpp(vts, vmarks, vbranching, *M0, *maxTime, *sims, *Bsamples, *Bfixed, *nuxiSteps, *CaSteps, *magSteps, *mu, *logC, *alpha, *nut, *xit, *numb, *ximb, *numa, *xima, vetasstepSds, vmagstepSds, vmuprior, *constrainOmori, vmus, vlogCs, valphas, vnuts, vxits, vnumbs, vximbs, vnumas, vximas,vBs);

		std::copy(vmus.begin(),vmus.end(),mus);
        std::copy(vlogCs.begin(),vlogCs.end(),logCs);
        std::copy(valphas.begin(),valphas.end(),alphas);
		std::copy(vnuts.begin(),vnuts.end(),nuts);
        std::copy(vxits.begin(),vxits.end(),xits);
		std::copy(vnumbs.begin(),vnumbs.end(),numbs);
        std::copy(vximbs.begin(),vximbs.end(),ximbs);
		std::copy(vnumas.begin(),vnumas.end(),numas);
        std::copy(vximas.begin(),vximas.end(),ximas);
		std::copy(vBs.begin(), vBs.end(), Bs);
		return;
	}

  // [[ export() ]]
  void estimateETAScorr_C(double *ts, double *marks, int *branching, int *n, double *M0, double *maxTime, int *sims, int *Bsamples, int *Bsize, int *Bfixed, int *nuxiSteps, int *CaSteps, int *magSteps,  double *mu, double *logC, double *alpha, double *nut, double *xit, double *numb, double *ximb, double *numa, double *xima, double *rho, double *etasstepSds, double *magstepSds, double *muprior, double *logCprior, int *constrainOmori, double *mus, double *logCs, double *alphas, double *nuts, double *xits, double *numbs, double *ximbs, double *numas, double *ximas, double *rhos, int *Bs){
		std::vector<double> vts(ts,ts + *n);
		std::vector<double> vmarks(marks,marks + *n);
		std::vector<int> vbranching(branching,branching + *n);
		std::vector<double> vetasstepSds(etasstepSds, etasstepSds + 4);
		std::vector<double> vmagstepSds(magstepSds, magstepSds + 5);
		std::vector<double> vmuprior(muprior, muprior + 2);
		std::vector<double> vlogCprior(logCprior, logCprior + 2);

   	 	std::vector<double> vmus;
    	std::vector<double> vlogCs;
    	std::vector<double> valphas;
   	 	std::vector<double> vnuts;
    	std::vector<double> vxits;
		std::vector<double> vnumbs;
		std::vector<double> vximbs;
		std::vector<double> vnumas;
		std::vector<double> vximas;
		std::vector<double> vrhos;
		std::vector<int> vBs;
    	vmus.reserve(*sims);
    	vlogCs.reserve(*sims);
    	valphas.reserve(*sims);
    	vnuts.reserve(*sims);
    	vxits.reserve(*sims);
		vnumbs.reserve(*sims);
    	vximbs.reserve(*sims);
		vnumas.reserve(*sims);
    	vximas.reserve(*sims);
    	vrhos.reserve(*sims);
		vBs.reserve(*Bsize);

		estimateETAScorr_Cpp(vts, vmarks, vbranching, *M0, *maxTime, *sims, *Bsamples, *Bfixed, *nuxiSteps, *CaSteps, *magSteps, *mu, *logC, *alpha, *nut, *xit, *numb, *ximb, *numa, *xima, *rho, vetasstepSds, vmagstepSds, vmuprior, vlogCprior, *constrainOmori, vmus, vlogCs, valphas, vnuts, vxits, vnumbs, vximbs, vnumas, vximas, vrhos, vBs);

		std::copy(vmus.begin(),vmus.end(),mus);
        std::copy(vlogCs.begin(),vlogCs.end(),logCs);
        std::copy(valphas.begin(),valphas.end(),alphas);
		std::copy(vnuts.begin(),vnuts.end(),nuts);
        std::copy(vxits.begin(),vxits.end(),xits);
		std::copy(vnumbs.begin(),vnumbs.end(),numbs);
        std::copy(vximbs.begin(),vximbs.end(),ximbs);
		std::copy(vnumas.begin(),vnumas.end(),numas);
        std::copy(vximas.begin(),vximas.end(),ximas);
        std::copy(vrhos.begin(),vrhos.end(),rhos);
		std::copy(vBs.begin(), vBs.end(), Bs);
		return;
	}

	// C functions to test Cpp functions against R
	// [[ export() ]]
	double dgpdlogC(double *x, double *nu, double *xi, double *mu){
		double value;
		value = exp(dgpdlog(*x, *nu, *xi, *mu));
		Rprintf("density: %f \n", value);
		return(value);
	}
  // [[ export() ]]
	double pgpdC(double *q, double *nu, double *xi, double *mu){
		double value = pgpd(*q, *nu, *xi, *mu);
		Rprintf("probability: %f \n", value);
		return(value);
	}
  // [[ export() ]]
	double gpd_to_gaussC(double *x_gpd, double *nu, double *xi, double *mu){
		double value = gpd_to_gauss(*x_gpd, *nu, *xi, *mu);
		Rprintf("x_gauss = %f \n", value);
		return(value);
	}
  // [[ export() ]]
	double gpdSurvivorC(double *x, double *nu, double *xi, double *mu){
		double value;
		value = gpdSurvivor(*x, *nu, *xi, *mu);
		Rprintf("probability: %f \n", value);
		return(value);
	}
  // [[ export() ]]
	void rnormC(double *samples, int *n, double *mu, double *sigma){
		std::vector<double> vsamples(samples, samples + *n);
		rnormCpp(vsamples, *n, *mu, *sigma);

		for (int i = 0; i < *n; i++){
		Rprintf("vsamples[%d] = %f \n", i, vsamples[i]);
		}

		std::copy(vsamples.begin(),vsamples.end(),samples);
		return;
	}
  // [[ export() ]]
	void discrete_sampleC(int *size, int *n, double *prob, int *samples){
		std::vector<int> vsamples(samples, samples + *n);
		std::vector<double> vprobs(prob,prob + *n);

		discrete_sampleCpp(*size, *n, vprobs, vsamples);

		for (int i = 0; i < *size; i++){
		Rprintf("vsamples[%d] = %d \n", i, vsamples[i]);
		}

		std::copy(vsamples.begin(),vsamples.end(),samples);
		return;
	}
  // [[ export() ]]
	void rdiscreteC(double *probs, int *n, int *sample){
		std::vector<double> vprobs(probs,probs + *n);
		int samp = -1;
		samp = rdiscrete(vprobs);
		*sample = samp;
		return;
	}
}


