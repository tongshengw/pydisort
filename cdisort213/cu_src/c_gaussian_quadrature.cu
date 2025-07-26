// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_gaussian_quadrature() =================*/

/*
   Compute weights and abscissae for ordinary Gaussian quadrature
   on the interval (0,1);  that is, such that
       sum(i=1 to M) ( GWT(i) f(GMU(i)) )
   is a good approximation to integral(0 to 1) ( f(x) dx )

   INPUT :     m        order of quadrature rule

   OUTPUT :    GMU(I)   array of abscissae (I = 1 TO M)
               GWT(I)   array of weights (I = 1 TO M)

   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
                 Integration, Academic Press, New York, pp. 87, 1975

   METHOD:     Compute the abscissae as roots of the Legendre polynomial P-sub-M using a cubically convergent
               refinement of Newton's method.  Compute the weights from eq. 2.7.3.8 of Davis/Rabinowitz.  Note
               that Newton's method can very easily diverge; only a very good initial guess can guarantee convergence.
               The initial guess used here has never led to divergence even for M up to 1000.

   ACCURACY:   Relative error no better than TOL or computer precision (DBL_EPSILON), whichever is larger

   INTERNAL VARIABLES:
    iter      : Number of Newton Method iterations
    pm2,pm1,p : 3 successive Legendre polynomials
    ppr       : Derivative of Legendre polynomial
    p2pri     : 2nd derivative of Legendre polynomial
    tol       : Convergence criterion for Legendre poly root iteration
    x,xi      : Successive iterates in cubically-convergent version of Newtons Method (seeking roots of Legendre poly.)

   Called by- c_dref, c_disort_set, c_surface_bidir
   Calls- c_errmsg
 -------------------------------------------------------------------*/

/* Maximum allowed iterations of Newton Method */
#define MAXIT 1000

void c_gaussian_quadrature(int    m,
                           double *gmu,
                           double *gwt)
{
  static int
    initialized = FALSE;
  register int
    iter,k,lim,nn,np1;
  double
    cona,t,en,nnp1,p=0,p2pri,pm1,pm2,ppr,
    prod,tmp,x,xi;
  static double
    tol;

  if (!initialized) {
    tol         = 10.*DBL_EPSILON;
    initialized = TRUE;
  }

  if (m < 1) {
    c_errmsg("gaussian_quadrature--Bad value of m",DS_ERROR);
  }

  if (m == 1) {
    GMU(1) = 0.5;
    GWT(1) = 1.0;
    return;
  }

  en   = (double)m;
  np1  = m+1;
  nnp1 = m*np1;
  cona = (double)(m-1)/(8*m*m*m);
  lim  = m/2;
  for (k = 1; k <= lim; k++) {
    /*
     * Initial guess for k-th root of Legendre polynomial, from Davis/Rabinowitz (2.7.3.3a)
     */
    t = (double)(4*k-1)*M_PI/(4*m+2);
    x = cos(t+cona/tan(t));

    /*
     * Upward recurrence for Legendre polynomials
     */
    for (iter = 1; iter <= MAXIT+1; iter++) {
      if (iter > MAXIT) {
        c_errmsg("gaussian_quadrature--max iteration count",DS_ERROR);
      }
      pm2 = 1.;
      pm1 = x;
      for (nn = 2; nn <= m; nn++) {
        p   = ((double)(2*nn-1)*x*pm1-(double)(nn-1)*pm2)/nn;
        pm2 = pm1;
        pm1 = p;
      }
      /*
       * Newton Method
       */
      tmp   = 1./(1.-x*x);
      ppr   = en*(pm2-x*p)*tmp;
      p2pri = (2.*x*ppr-nnp1*p)*tmp;
      xi    = x-p/ppr*(1.+p/ppr*p2pri/(2.*ppr));
      /*
       * Check for convergence
       */
      if (fabs(xi-x) <= tol) {
        break;
      }
      else {
        x = xi;
      }
    }

    /*
     * Iteration finished--calculate weights, abscissae for (-1,1)
     */
    GMU(k)     = -x;
    GWT(k)     = 2./(tmp*SQR(en*pm2));
    GMU(np1-k) = -GMU(k);
    GWT(np1-k) =  GWT(k);
  }

  /*
   * Set middle abscissa and weight for rules of odd order
   */
  if (m%2 != 0) {
    GMU(lim+1) = 0.;
    prod       = 1.;
    for (k = 3; k <= m; k+=2) {
      prod *= (double)k/(k-1);
    }
    GWT(lim+1) = 2./SQR(prod);
  }
  /*
   * Convert from (-1,1) to (0,1)
   */
  for (k = 1; k <= m; k++) {
    GMU(k) = 0.5*GMU(k)+0.5;
    GWT(k) = 0.5*GWT(k);
  }

  return;
}

#undef MAXIT

/*============================= end of c_gaussian_quadrature() ==========*/
