// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_legendre_poly() =======================*/

/*
       Computes the normalized associated Legendre polynomial, defined
       in terms of the associated Legendre polynomial Plm = P-sub-l-super-m as

          Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

       for fixed order m and all degrees from l = m to TWONM1.
       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
       from a prior call to the routine.

       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of High-Order
                    Associated Legendre Polynomials, J. Quant. Spectrosc. Radiat. Transfer 10,
                    557-562, 1970. (hereafter D/A)

       METHOD: Varying degree recurrence relationship.

       NOTES:
       (1) The D/A formulas are transformed by setting m=n-1; l=k-1.
       (2) Assumes that routine is called first with  m = 0, then with
           m = 1, etc. up to  m = twonm1.


  I N P U T     V A R I A B L E S:

       nmu    :  Number of arguments of YLM
       m      :  Order of YLM
       maxmu  :
       twonm1 :  Max degree of YLM
       MU(i)  :  Arguments of YLM (i = 1 to nmu)

       If m > 0, YLM(m-1,i) for i = 1 to nmu is assumed to exist from a prior call.


  O U T P U T     V A R I A B L E:

       YLM(l,i) :  l = m to twonm1, normalized associated Legendre polynomials
                   evaluated at argument MU(i)

   Called by- c_disort, c_albtrans
 -------------------------------------------------------------------*/

void c_legendre_poly(int     nmu,
                     int     m,
                     int     maxmu,
                     int     twonm1,
                     double *mu,
                     double *ylm)
{
  register int
    i,l;
  register double
    tmp1,tmp2;

  if (m == 0) {
    /*
     * Upward recurrence for ordinary Legendre polynomials
     */
    for (i = 1; i <= nmu; i++) {
      YLM(0,i) = 1.;
      YLM(1,i) = MU(i);
    }
    for (l = 2; l <= twonm1; l++) {
      for (i = 1; i <= nmu; i++) {
        YLM(l,i) = ((double)(2*l-1)*MU(i)*YLM(l-1,i)-(double)(l-1)*YLM(l-2,i))/l;
      }
    }
  }
  else {
    for (i = 1; i <= nmu; i++) {
      /*
       * Y-sub-m-super-m; derived from D/A eqs. (11,12), STWL(58c)
       */
      YLM(m,i) = -sqrt((1.-1./(2*m))*(1.-SQR(MU(i))))*YLM(m-1,i);

      /*
       * Y-sub-(m+1)-super-m; derived from D/A eqs.(13,14) using eqs.(11,12), STWL(58f)
       */
      YLM(m+1,i) = sqrt(2.*m+1.)*MU(i)*YLM(m,i);
    }
    /*
     * Upward recurrence; D/A eq.(10), STWL(58a)
     */
    for (l = m+2; l <= twonm1; l++) {
      tmp1 = sqrt((l-m  )*(l+m  ));
      tmp2 = sqrt((l-m-1)*(l+m-1));
      for (i = 1; i <= nmu; i++) {
        YLM(l,i) = ((double)(2*l-1)*MU(i)*YLM(l-1,i)-tmp2*YLM(l-2,i))/tmp1;
      }
    }
  }

  return;
}

/*============================= end of c_legendre_poly() ================*/
