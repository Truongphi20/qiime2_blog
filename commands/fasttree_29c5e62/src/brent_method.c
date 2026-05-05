#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/******************************************************************************/
/* Minimization of a 1-dimensional function by Brent's method (Numerical Recipes)            
 * Borrowed from Tree-Puzzle 5.1 util.c under GPL
 * Modified by M.N.P to pass in the accessory data for the optimization function,
 * to use 2x bounds around the starting guess and expand them if necessary,
 * and to use both a fractional and an absolute tolerance
 */

#define ITMAX 100
#define CGOLD 0.3819660
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent(double ax, double bx, double cx, double (*f)(double, void *), void *data,
	     double ftol, double atol,
	     double *foptx, double *f2optx, double fax, double fbx, double fcx)
{
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if (fax < fcx) {
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else {
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;	
	}
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=ftol*fabs(x);
		tol2=2.0*(tol1+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))
		    || fabs(a-b) < atol) {
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
				(v*v*xw + x*x*wv + w*w*vx);
			return x;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,data);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
		(v*v*xw + x*x*wv + w*w*vx);
	return x;
} /* brent */
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN

/* one-dimensional minimization - as input a lower and an upper limit and a trial
  value for the minimum is needed: xmin < xguess < xmax
  the function and a fractional tolerance has to be specified
  onedimenmin returns the optimal x value and the value of the function
  and its second derivative at this point
  */
double onedimenmin(double xmin, double xguess, double xmax, double (*f)(double,void*), void *data,
		   double ftol, double atol,
		   /*OUT*/double *fx, /*OUT*/double *f2x)
{
	double optx, ax, bx, cx, fa, fb, fc;
		
	/* first attempt to bracketize minimum */
	if (xguess == xmin) {
	  ax = xmin;
	  bx = 2.0*xguess;
	  cx = 10.0*xguess;
	} else if (xguess <= 2.0 * xmin) {
	  ax = xmin;
	  bx = xguess;
	  cx = 5.0*xguess;
	} else {
	  ax = 0.5*xguess;
	  bx = xguess;
	  cx = 2.0*xguess;
	}
	if (cx > xmax)
	  cx = xmax;
	if (bx >= cx)
	  bx = 0.5*(ax+cx);
	if (verbose > 4)
	  fprintf(stderr, "onedimenmin lo %.4f guess %.4f hi %.4f range %.4f %.4f\n",
		  ax, bx, cx, xmin, xmax);
	/* ideally this range includes the true minimum, i.e.,
	   fb < fa and fb < fc
	   if not, we gradually expand the boundaries until it does,
	   or we near the boundary of the allowed range and use that
	*/
	fa = (*f)(ax,data);
	fb = (*f)(bx,data);
	fc = (*f)(cx,data);
	while(fa < fb && ax > xmin) {
	  ax = (ax+xmin)/2.0;
	  if (ax < 2.0*xmin)	/* give up on shrinking the region */
	    ax = xmin;
	  fa = (*f)(ax,data);
	}
	while(fc < fb && cx < xmax) {
	  cx = (cx+xmax)/2.0;
	  if (cx > xmax * 0.95)
	    cx = xmax;
	  fc = (*f)(cx,data);
	}
	optx = brent(ax, bx, cx, f, data, ftol, atol, fx, f2x, fa, fb, fc);

	if (verbose > 4)
	  fprintf(stderr, "onedimenmin reaches optimum f(%.4f) = %.4f f2x %.4f\n", optx, *fx, *f2x);
	return optx; /* return optimal x */
} /* onedimenmin */