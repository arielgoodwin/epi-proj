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
const int arrsize = 1000;

// function declarations
double gradTheta(double lambda, double x[],int* g);
double newton(double y[], double stopThr, double guess);

int main() {

		// code for testing the algorithm numerically //

		// initialize variables
		double *y;
		double *x;
		srand(time(NULL));
		int 	i,j;
	    int 	Nbrea;
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
		double lambda_warm = -1;

		// for each test example
		for (j=0; j<=Nbrea; j++) {

			// initialize the vector y with random values from a given probability distribution
			// here we use gaussian with sigma = 0.1
			for (i=0; i<arrsize; i++) {
		    	y[i]=gsl_ran_gaussian(r,0.1);
		    }

			// choose a random entry and increase it artificially to ensure ||y||_1 > 1
		    y[(int)(rand() / (((double)RAND_MAX+1.0)/arrsize))]+=a;

		    // time the algorithm
		    start=clock();

		    // algorithm function call
		    // here we have implemented the warm-start version, the function can be modified accordingly to use a different initialization
		    // if desired
		    lambda_warm = newton(y,0.001,lambda_warm);

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

// this is the derivative (w.r.t. lambda) of the Theta function in Hoheisel's paper
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
double newton(double *y, double stopThr, double guess) {

	// need an initial point to start the method
	// we use the following method: sample sqrt(N)*log(N) random coordinates
	// of the vector we are projecting, and take the largest of their absolute values
	// as the initial guess. in the warm start version, this method is only used in the first trial.

	double lambda_k;
	if (guess < 0){
		double mx = 0;

		for (int i = 0; i < log(arrsize)*sqrt(arrsize); i++){
			double val = 0;
			if ((val = fabs(y[(int)(rand() / (((double)RAND_MAX+1.0)/arrsize))])) > mx){
				mx = val;
			}
		}

		lambda_k = mx;

	} else{
		lambda_k = guess;
	}


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

	return lambda_k;

}
