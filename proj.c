/*
 * proj.c
 *
 * This program uses a semismooth Newton method to
 * compute the projection of vectors in n-dimensional space onto
 * the l1 ball.
 *
 *  Created on: May 28, 2020
 *      Author: ariel goodwin
 */

// useful includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// dimension we are working in
const int arrsize = 20;

// function declarations
double Theta(double lambda, double x[]);
double gradTheta(double lambda, double x[],int* g);
void projball(double y[], double stopThr);
void newton(double y[], double stopThr);

int main() {

		// code for testing the algorithm numerically //

		// initialize variables
		double *y;
		double *x;
		srand(time(NULL));
		int 	i,j;
	    int 	Nbrea;
	   // unsigned int length;
	    clock_t	start, end;
	    const gsl_rng_type *T;
	  	gsl_rng *r;
	  	gsl_rng_env_setup();
		gsl_rng_default_seed=rand();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		double timemean=0.0, timevar=0.0, timedelta;
		const double a = 1.0;
		srand((unsigned int)pow(time(NULL)%100,3));
	    double*	 timetable;

	    // a vector that can be used to store the projected vector if its undesirable to overwrite the original vector
	    x=(double*)malloc(arrsize*sizeof(double));

	    // how many test examples you want to run
		Nbrea = 10000;


		y=(double*)malloc(arrsize*sizeof(double));
		timetable=(double*)malloc((Nbrea+1)*sizeof(double));

		// for each test example
		for (j=0; j<=Nbrea; j++) {

			// initialize the vector y with random values from a given probability distribution
			// here we use gaussian with sigma = 0.1
			for (i=0; i<arrsize; i++) {
		    	y[i]=gsl_ran_gaussian(r,0.1);
		    }

			// choose a random entry and increase it artificially to ensure ||y||_1 > 1
		   // y[(int)(rand() / (((double)RAND_MAX+1.0)/arrsize))]+=a;

		    // time the algorithm
		    start=clock();

		    // algorithm function call
		    newton(y,0.001);

		    end=clock();
		    timetable[j]=(double)(end-start)/CLOCKS_PER_SEC;
		}

		// we discard the first value, because the computation time is higher,
		//probably because of cache operations
		for (j=1; j<=Nbrea; j++) {
			timedelta=timetable[j]-timemean;
			timemean+=timedelta/j;
			timevar+=timedelta*(timetable[j]-timemean);
		}
		timevar=sqrt(timevar/Nbrea);

		// print the average time
		printf("av. time: %e std dev: %e\n",timemean,timevar);
	    return 0;

}

// this is the derivative (w.r.t. lambda) of the Theta function below
// we do some extra work in this function to compute a Bouligand subgradient at lambda_k
// the code is written assuming the ball has unit radius
double gradTheta(double lambda, double x[], int* g) {

	// variables to keep track of computed values
	double sum = 0;
	int count = 0;
	double ab = 0;
	int i;

	// iterate through the vector
	for (i = 0; i < arrsize; i++) {

		// add thresholded entry to the sum
		ab = fabs(x[i]);
		sum += fmax(ab - lambda, 0);

		// check if this entry contributes to the Bouligand subdifferential
		if (ab > lambda) {
			count++;
		}
	}

	// set this pointer to the value of the Bouligand subgradient, will be used in the larger function
	*g = count;

	// return the final value or 0 if the optimal value is at the boundary (lambda = 0)
	if (lambda == 0 && 1-sum > 0){
		return 0;
	} else {
		return 1 - sum;
	}

}

// newton's method - a slight variation of [Hoheisel, Algorithm 1]
// this algorithm does not use a linesearch - it can be shown that once the iterates are smaller than lambda_opt,
// the algorithm converges monotonically (increasing) to the root.
void newton(double *y, double stopThr) {

	// need an initial point to start the method
	// we use the following method: sample sqrt(N)*log(N) random coordinates
	// of the vector we are projecting, and take the largest of their absolute values
	// as the initial guess.

	double mx = 0;

	for (int i = 0; i < log(arrsize)*sqrt(arrsize); i++){
		double val = 0;
		if ((val = fabs(y[(int)(rand() / (((double)RAND_MAX+1.0)/arrsize))])) > mx){
			mx = val;
		}
	}

	double lambda_k = mx;

	// these variables hold information so that less needs to be recomputed
	double F;
	int boul = 0;

	// iterate while checking termination condition
	// note that a value of the Bouligand subdifferential is computed in this step
	while (fabs(F=gradTheta(lambda_k,y,&boul)) >= stopThr ){

		// take a full newton step and reset the subdifferential
		lambda_k += fmax(-lambda_k, -F/boul);
		boul = 0;

	}

	// compute projection using optimal lambda
	for (int i=0; i<arrsize; i++){
		y[i]=(y[i]-lambda_k>0.0 ? y[i]-lambda_k : (y[i]+lambda_k<0.0 ? y[i]+lambda_k : 0.0));
	}

}

// the Theta function defined in [Hoheisel, Corollary 4.6]
// the code is written assuming the l1 ball has unit radius
double Theta(double lambda, double x[]) {

	// variables for holding sums
	double sum = 0;
	double sumsqr = 0;

	// iterate through all entries to compute the sums
	for (int i = 0; i < arrsize; i++) {

		sum += fmax(fabs(x[i])-lambda,0);
		sumsqr += (x[i]-copysign(1.0,x[i])*fmax(fabs(x[i]) - lambda,0))*((x[i]-copysign(1.0,x[i])*fmax(fabs(x[i]) - lambda,0)));

	}

	// final computation: lamb + phibar_{x,f}(lamb)
	return lambda*(1-sum) - 0.5*sumsqr;

}

// this function is more faithful to [Hoheisel, Algorithm 1] - it actually uses a linesearch
// to ensure global convergence. however, based on the convergence results we've proved for
// the first newton method, it is only included here to give an idea of how Algorithm 1 may be
// implemented for other applications
void projball(double *y, double stopThr) {

	// termination condition
	double delta = stopThr;

	// armijo search parameters
	double sigma = 0.0001;
	double beta = 0.5;

	// random coordinate initial guess method
	double mx = 0;

	for (int i = 0; i < log(arrsize)*sqrt(arrsize); i++){
		double val = 0;
		if ((val = fabs(y[(int)(rand() / (((double)RAND_MAX+1.0)/arrsize))])) > mx){
			mx = val;
		}
	}

	// initial guess
	double lambda_k = mx;

	// these variables hold information so that less needs to be recomputed
	double F;
	double d_k;
	double p = 1;
	int boul = 0;

	// iterate while checking termination condition
	while (fabs(F=gradTheta(lambda_k,y,&boul)) >= delta){

		// used to control the stepsize
		p = 1;

		// find the descent direction
		d_k = fmax(-lambda_k, -F/boul);

		// armijo line search
		while (Theta(lambda_k + p*d_k,y) > Theta(lambda_k,y) + p*sigma*F*d_k) {

			p*=beta;

		}

		// update test point based on line search stepsize
		lambda_k = lambda_k + p*d_k;
	}

	// compute projection using optimal lambda
	for (int i=0; i<arrsize; i++){
		y[i]=(y[i]-lambda_k>0.0 ? y[i]-lambda_k : (y[i]+lambda_k<0.0 ? y[i]+lambda_k : 0.0));
	}

}




