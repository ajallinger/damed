/* test program */
#include <stdio.h>
#include <math.h>
int damed (double *y, double *dydx, int n, double tau,
            double *ymed, double *dmdx);
#define N 10  /* problem size */

int main (void)
{
/* Initialized data */
	static double x[N] = { 3.,8.,3.,-1.,0.,5.,4.,8.,9.,2. };
	static double dydx[N] = { 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. };
/* Local variables */
	int ierr, j;
	double d, y[N], ave, dev, tau, dmdx, ymed;
/* Begin. */
	ave = 0.;
	for (j = 0; j < N; ++j) {
		y[j] = x[j];
		ave += y[j];
	}
	ave /= (double) N;
	dev = 0.;
	for (j = 0; j < N; ++j) {
		d = y[j] - ave;
		dev += d * d;
	}
	dev = sqrt(dev / (double) N);
	printf("x is: ");
	for (j = 0; j < N; ++j) printf(" %f", x[j]);
	printf("\n");
	tau = dev / 10.;
	printf("using tau: %f\n", tau);
	j = (int) N;
	ierr = damed (y, dydx, j, tau, &ymed, &dmdx);
	printf("damed returns #%i median: %f\n", ierr, ymed);
	printf("derivative: %f\n", dmdx);
	return 0;
} /* end of main */
