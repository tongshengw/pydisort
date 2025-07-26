// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_planck_func2() =========================*/

/*
  Computes Planck function integrated between two wavenumbers,
  except if wnmulo = wnmuhi, then the Planck function at wnumlo is returned

  I N P U T :  wnumlo : Lower wavenumber [inv cm] of spectral interval
               wnumhi : Upper wavenumber
               t      : Temperature [K]

  O U T P U T :  ans  : Integrated Planck function [Watts/sq m]
                         = integral (wnumlo to wnumhi) of 2h c*c nu*nu*nu/(exp(hc nu/(kT))-1),
                         where h = Plancks constant, c = speed of light, nu = wavenumber,
                         T=temperature,and k = Boltzmann constant

  REFERENCE : Specifications of the physical world: New value of the fundamental constants,
                Dimensions/N.B.S., Jan. 1974

  METHOD :  For  -wnumlo-  close to  -wnumhi-, a Simpson-rule quadrature is done
            to avoid ill-conditioning; otherwise

            (1)  For wavenumber (wnumlo or wnumhi) small, integral(0 to wnum) is calculated by expanding
                 the integrand in a power series and integrating term by term;

            (2)  Otherwise, integral(wnumlo/hi to infinity) is calculated by expanding the denominator of the
                 integrand in powers of the exponential and integrating term by term.

  ACCURACY :  At least 6 significant digits, assuming the physical constants are infinitely accurate

  ERRORS that are not trapped:

      * Power or exponential series may underflow, giving no significant digits.
        This may or may not be of concern, depending on the application.

      * Simpson-rule special case is skipped when denominator of integrand will cause overflow.
        In that case the normal procedure is used, which may be inaccurate if the wavenumber limits
        (wnumlo, wnumhi) are close together.
 ----------------------------------------------------------------------

        LOCAL VARIABLES

        a1,2,... :  Power series coefficients
        c2       :  h*c/k, in units cm*k (h = Planck's constant, c = speed of light, k = Boltzmann constant)
        D(I)     :  Exponential series expansion of integral of Planck function from wnumlo (i=1)
                    or wnumhi (i=2) to infinity
        ex       :  exp(-V(I))
        exm      :  pow(ex,m)
        mmax     :  No. of terms to take in exponential series
        mv       :  multiples of 'V(i)'
        P(I)     :  Power series expansion of integral of Planck function from zero to wnumlo (i=1) or wnumhi (i=2)
        sigma    :  Stefan-Boltzmann constant (W m-2 K-4)
        sigdpi   :  sigma/pi
        smallv   :  Number of times the power series is used (0,1,2)
        V(I)     :  c2*(wnumlo(i=1) or wnumhi(i=2))/temperature
        vcut     :  Power-series cutoff point
        vcp      :  Exponential series cutoff points
        vmax     :  Largest allowable argument of 'exp' function
  ----------------------------------------------------------------------*/

#define A1    (1./3.)
#define A2    (-1./8.)
#define A3    (1./60.)
#define A4    (-1./5040.)
#define A5    (1./272160.)
#define A6    (-1./13305600.)
#define C2    (1.438786)
#define SIGMA (5.67032e-8)
#define VCUT  (1.5)
#define PLKF(x) ({const double _x = (x); _x*_x*_x/(exp(_x)-1.);})

double __attribute__((weak)) c_planck_func2(double wnumlo,
                      double wnumhi,
                      double t)
{
  register int
    m,n,smallv,k,i,mmax;
  static int
    initialized = FALSE;
  double
    ans,del,val,val0,oldval,exm,
    ex,mv,vsq,wvn,arg,hh,
    d[2],p[2],v[2];
  const double
    vcp[7] = {10.25,5.7,3.9,2.9,2.3,1.9,0.0};
  static double
    sigdpi,vmax,conc,c1;

  if (!initialized) {
    sigdpi = SIGMA/M_PI;
    vmax   = log(DBL_MAX);
    conc   = 15./pow(M_PI,4.);
    c1     = 1.1911e-8;

    initialized = TRUE;
  }
  if (t < 0. || wnumhi < wnumlo || wnumlo < 0.) {
    c_errmsg("planck_func2--temperature or wavenumbers wrong",DS_ERROR);
  }
  if (t < 1.e-4) {
    return 0.;
  }
  if (wnumhi == wnumlo) {
    wvn    = wnumhi;
    arg    = exp(-C2*wvn/t);
    return c1*wvn*wvn*wvn*arg/(1.-arg);
  }

  v[0] = C2*wnumlo/t;
  v[1] = C2*wnumhi/t;

  if (v[0] > DBL_EPSILON && v[1] < vmax && (wnumhi-wnumlo)/wnumhi < 1.e-2) {
    /*
     * Wavenumbers are very close. Get integral by iterating Simpson rule to convergence.
     */
    hh     = v[1]-v[0];
    oldval = 0.;
    val0   = PLKF(v[0])+PLKF(v[1]);
    for (n = 1; n <= 10; n++) {
      del = hh/(2*n);
      val = val0;
      for (k = 1; k <=2*n-1; k++) {
        val += (double)(2*(1+k%2))*PLKF(v[0]+(double)k*del);
      }
      val *= del*A1;
      if (fabs((val-oldval)/val) <= 1.e-6) {
        return sigdpi*SQR(t*t)*conc*val;
      }
      oldval = val;
    }
    c_errmsg("planck_func2--Simpson rule did not converge",DS_WARNING);
    return sigdpi*SQR(t*t)*conc*val;
  }

  smallv = 0;
  for (i = 0; i <= 1; i++) {
    if(v[i] < VCUT) {
      /*
       * Use power series
       */
      smallv++;
      vsq  = v[i]*v[i];
      p[i] = conc*vsq*v[i]*(A1+v[i]*(A2+v[i]*(A3+vsq*(A4+vsq*(A5+vsq*A6)))));
    }
    else {
      /*
       * Use exponential series
       *
       * Find upper limit of series
       */
      mmax = 1;
      while (v[i] < vcp[mmax-1]) {
        mmax++;
      }

      ex   = exp(-v[i]);
      exm  = 1.;
      d[i] = 0.;

      for (m = 1; m <= mmax; m++) {
        mv    = (double)m*v[i];
        exm  *= ex;
        d[i] += exm*(6.+mv*(6.+mv*(3.+mv)))/SQR(m*m);
      }
      d[i] *= conc;
    }
  }

  if (smallv == 2) {
    /*
     * wnumlo and wnumhi both small
     */
    ans = p[1]-p[0];
  }
  else if (smallv == 1) {
    /*
     * wnumlo small, wnumhi large
     */
    ans = 1.-p[0]-d[1];
  }
  else {
    /*
     * wnumlo and wnumhi both large
     */
    ans = d[0]-d[1];
  }
  ans *= sigdpi*SQR(t*t);
  if (ans == 0.) {
    c_errmsg("planck_func2--returns zero; possible underflow",DS_WARNING);
  }

  return ans;
}

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef C2
#undef SIGMA
#undef VCUT
#undef PLKF

/*============================= end of c_planck_func2() =================*/
