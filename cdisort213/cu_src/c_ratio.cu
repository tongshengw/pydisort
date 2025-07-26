// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_ratio() ===============================*/

/*
 * Calculate ratio a/b with overflow and underflow protection
 * (thanks to Prof. Jeff Dozier for some suggestions here).
 *
 * Modification in this C version: in the case b == 0., returns 1.+a.
 *
 * Called by: c_disort
 */

double c_ratio(double a,
             double b)
{
  static int
    initialized = FALSE;
  static double
    tiny,huge,powmax,powmin;
  double
    ans,absa,absb,powa,powb;

  if(!initialized) {
    tiny   = DBL_MIN;
    huge   = DBL_MAX;
    powmax = log10(huge);
    powmin = log10(tiny);

    initialized = TRUE;
  }

  if (c_fcmp(b,0.) == 0) {
    ans = 1.+a;
  }
  else if (c_fcmp(a,0.) == 0) {
    ans = 0.;
  }
  else {
    absa = fabs(a);
    absb = fabs(b);
    powa = log10(absa);
    powb = log10(absb);
    if (c_fcmp(absa,tiny) < 0 && c_fcmp(absb,tiny) < 0) {
      ans = 1.;
    }
    else if (c_fcmp(powa-powb,powmax) >= 0) {
      ans = huge;
    }
    else if(c_fcmp(powa-powb,powmin) <= 0) {
      ans = tiny;
    }
    else {
      ans = absa/absb;
    }

   /*
    * NOTE: Don't use old trick of determining sign from a*b because a*b
    *       may overflow or underflow.
    */
    if ( (a > 0. && b < 0.) || (a < 0. && b > 0.) ) {
      ans *= -1;
    }
  }

  return ans;
}

/*============================= end of c_ratio() ========================*/
