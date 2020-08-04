/*
 * proj.c
 *
 * This program uses a semismooth Newton method to
 * compute the projection of vectors in n-dimensional space onto
 * the epigraph of negative sum of logs function.
 *
 *  Created on: May 28, 2020
 *      Author: arielgoodwin
 */

// useful includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// dimension we are working in
const int arrsize = 1;

// function declarations
double gradTheta(double lambda, double x[], double alpha);
double gradgradTheta(double lambda, double x[]);
void projlognewt(double *y, double alpha, double stopThr);


int main() {

	// code for testing the functions numerically
	srand(time(NULL));
	double*	y;
	double*	x;
	double alpha;
	int 	i,j;
	int 	Nbrea;
    unsigned int length;
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

    length = arrsize;

    // number of trials you wish to run
    Nbrea = 10000;
	y=(double*)malloc(length*sizeof(double));

	// a second vector if it is desireable to compute the projection without erasing original vector
	x=(double*)malloc(length*sizeof(double));
	timetable=(double*)malloc((Nbrea+1)*sizeof(double));

	// for every trial
	for (j=0; j<=Nbrea; j++) {

		// initialize a random vector with coordinates in uniform distribution over [-1,1]
		for (i=0; i<length; i++) {

			y[i]=gsl_ran_flat(r,-1,1);

	    }

		// initialize random real value alpha uniformly over [-2,-0.5] (this is the second coordinate of the point being projected)
		alpha = gsl_ran_flat(r,-2,-0.5);
	    start=clock();

	    // function call
	    projlognewt(y,alpha,0.001);

	    end=clock();
	    timetable[j]=(double)(end-start)/CLOCKS_PER_SEC;
	}
	// we discard the first value, because the computation time is higher,
	//	probably because of cache operations
	for (j=1; j<=Nbrea; j++) {
		timedelta=timetable[j]-timemean;
		timemean+=timedelta/j;
		timevar+=timedelta*(timetable[j]-timemean);
	}
	timevar=sqrt(timevar/Nbrea);
	printf("av. time: %e std dev: %e\n",timemean,timevar);
	free(x);
	free(y);
    return 0;

}

// derivative of Theta with w.r.t. lambda (Theta being the function from Hoheisel's paper)
double gradTheta(double lambda, double x[], double alpha) {

	double sum = 0;

	for (int i = 0; i < arrsize; i++) {

			sum += log((x[i]+sqrt(x[i]*x[i]+4*lambda))/2.0);

	}

	return lambda + alpha + sum;

}

// the second derivative of Theta w.r.t. lambda
double gradgradTheta(double lambda, double x[]) {

	double sum = 0;

	for (int i = 0; i < arrsize; i++) {

		sum += 1.0/(x[i]*sqrt(x[i]*x[i]+4*lambda)+x[i]*x[i]+4*lambda);

	}

	return 1 + 2*sum;

}

// function that does the actual projection
void projlognewt(double *y, double alpha, double stopThr) {

		// termination condition
		double delta = stopThr;

		// initial guess and iteration counter (could employ a warm start here too)
		double lambda_k = sqrt(arrsize);
		int k = 0;

		// these variables hold information so that less needs to be recomputed
		double F;
		double g_k;
		double d_k;
		double epsilon_k;

		// iterate while checking termination condition
		while ((fabs(F=gradTheta(lambda_k,y,alpha)) >= delta)) {

			// small positive sequence going to zero
			epsilon_k = 1.0/pow(2.0,k+1);

			// set the descent direction appropriately, keeping away from 0
			if (F < 0) {

				g_k = gradgradTheta(lambda_k,y);
				d_k = -F/g_k;


			} else {

				g_k = gradgradTheta(lambda_k,y);
				d_k = fmax(-lambda_k + epsilon_k, -F/g_k);


			}

			// update test point and increment index
			lambda_k = lambda_k + d_k;
			k++;
		}

		// compute projection
		for (int i = 0; i < arrsize; i++){

			y[i] = (y[i] + sqrt(y[i]*y[i]+4*lambda_k))/2.0;

		}

}
