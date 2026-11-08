/*----------------------------------------------------------------------
disect - solve an equation by bisection, allowing an interval f = 0.
by Andy Allinger, 2018-2023, released to the public domain.

   Permission  to  use, copy, modify, and distribute this software and
   its documentation  for  any  purpose  and  without  fee  is  hereby
   granted,  without any conditions or restrictions.  This software is
   provided "as is" without express or implied warranty.

___Name______Type_______In/Out___Description____________________________
   func      function   In       The function to be solved
                                  It has the calling sequence:
                                   f = func (x, darg, n, tau)
   b         double*    Both     One end of the interval to search
   c         double*    Both     The other end of the interval
                                  On return, B and C have been changed
                                  to bounds on the solution.
   atol      double     In       Absolute error tolerance for solution
   darg[n]   double*    Thru     Array passed to user function
   n         int        Thru     Length of darg
   tau       double     Thru     Argument passed to user function

   Return value:                 0 = no errors
                                 1 = no zero crossing
                                 2 = too many iterations
----------------------------------------------------------------------*/
int disect (double (*func)(double, double*, int, double),
        double *b, double *c, double atol, double *darg, int n, double tau)
{
	const int MAXIT = 9000;
	int kount;
	double f, fb, fc, xg, xl, rootg, rootl;
/* Function Body */
	if (*b > *c) {  /* Ensure that b < c */
		f = *b;
		*b = *c;
		*c = f;
	}
	if (*c - *b < atol) return 0;
/* Check that a root is present */
	fb = (*func)(*b, darg, n, tau);
	fc = (*func)(*c, darg, n, tau);
	kount = 2;
if (! ((fb < 0. && fc > 0.) || (fb > 0. && fc < 0.)) ) return 1;  /* No change of sign */
/* Find rootl, the least value in [b,c] such that f(rootl) = 0 */
	xl = *b;
	xg = *c;
L10:
	rootl = 0.5 * (xl + xg);
	if (xg - xl <= atol) goto L20;
	f = (*func)(rootl, darg, n, tau);
	++kount;
	if (kount > MAXIT) return 2;
	if (fb < 0.) {  /* search for f(x) < 0 */
		if (f < 0.) {
			xl = rootl;
		} else {
			xg = rootl;
		}
	} else {  /* search for f(x) > 0 */
		if (f > 0.) {
			xl = rootl;
		} else {
			xg = rootl;
		}
	}
	goto L10;
/* Find rootg, the greatest value in [b,c] such that f(rootg) = 0 */
L20:
	xl = *b;
	xg = *c;
L30:
	rootg = 0.5 * (xl + xg);
	if (xg - xl <= atol) goto L40;
	f = (*func)(rootg, darg, n, tau);
	++kount;
	if (kount > MAXIT) return 2;
	if (fc < 0.) {  /* search for f(x) < 0 */
		if (f < 0.) {
			xg = rootg;
		} else {
			xl = rootg;
		}
	} else {  /* search for f(x) > 0 */
		if (f > 0.) {
			xg = rootg;
		} else {
			xl = rootg;
		}
	}
	goto L30;
/* Final answer */
L40:
	*b = rootl;
	*c = rootg;
	return 0;
} /* end of disect */
