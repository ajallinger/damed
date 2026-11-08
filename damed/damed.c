/*----------------------------------------------------------------------
damed.c - differentiable approximation to the median
by Andy Allinger, 2023, released to the public domain

   Permission  to  use, copy, modify, and distribute this software and
   its documentation  for  any  purpose  and  without  fee  is  hereby
   granted,  without any conditions or restrictions.  This software is
   provided "as is" without express or implied warranty.

Intended for use in unconstrained optimization.
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <float.h>
#include <math.h>

bool dafdiv (double numer, double denom);
int disect (double (*func)(double, double*, int, double),
        double *b, double *c, double atol, double *darg, int n, double tau);

/*----------------------------------------------------------------------
              The counting function, using approximate absolute value
----------------------------------------------------------------------*/
double countf (double ymed, double *y, int n, double tau)
{
/* Local variables */
	int j;
	double c, d, oot, s, t, xj;
/* Function Body */
	oot = 1. / tau;
	c = 0.;  /* compensated addition */
	s = 0.;
	for (j = 0; j < n; ++j) {
		d = (y[j] - ymed) * oot;
		xj = erf(d);
		t = s + xj;
		if (fabs(s) >= fabs(xj)) {
			c += (s - t) + xj;
		} else {
			c += (xj - t) + s;
		}
		s = t;
	}
	d = s + c;
	return d;
} /* end of countf */

/*----------------------------------------------------------------------
damed - differentiable approximation to the median

___Name______Type_______In/Out___Description____________________________
   y[n]      double*    In       Array of functions of x
   dydx[n]   double*    In       Partial derivatives dy[i]/dx
   n         int        In       Length of Y
   tau       double     In       Small number
   ymed      double*    Out      Smooth approximate median of Y
   dmdx      double*    Out      Derivative of the smooth median

   Return value                  0 = no errors
                                 1 = division by zero
                                 2 = errors from root finding
----------------------------------------------------------------------*/
int damed (double *y, double *dydx, int n, double tau,
            double *ymed, double *dmdx)
{
/* Local variables */
	int iflag, j;
	double atol, b, bel, c, d, oot, r, s, z, zmin;
/* Function Body */
	b = y[0];  /* range of Y */
	c = y[0];
	for (j = 1; j < n; ++j) {
		if (y[j] < b) b = y[j];
		if (y[j] > c) c = y[j];
	}
	r = (c + b) * 0.5;
	s = (c - b) * 0.5;
	atol = DBL_EPSILON * fmax(r,s);
	iflag = disect (countf, &b, &c, atol, y, n, tau);
	if (iflag != 0) return 2;
	*ymed = 0.5 * (b + c);
	oot = 1. / tau;
	zmin = 1e36;
	for (j = 0; j < n; ++j) {
		d = (y[j] - *ymed) * oot;
		z = d * d;
		if (z < zmin) zmin = z;
	}
	r = 0.;
	s = 0.;
	for (j = 0; j < n; ++j) {
		d = (y[j] - *ymed) * oot;
		bel = exp(zmin - d * d);
		r += bel * dydx[j];
		s += bel;
	}
	*dmdx = 0.;
	if (!dafdiv(r, s)) return 1;
	*dmdx = r / s;
	return 0;
} /* end of damed */
