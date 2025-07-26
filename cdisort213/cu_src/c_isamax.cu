// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_isamax() ===============================*/

/*
 INPUT--  n        Number of elements in vector of interest
          sx       Array, length n, containing vector

 OUTPUT-- ans      First i, i = 1 to n, to maximize fabs(SX(i))

 NOTE: Fortran input incx removed because it is not used by
       disort or twostr
 ---------------------------------------------------------------------*/

int c_isamax(int     n,
             double *sx)
{
  register int
    ans=0,i;
  double
   smax,xmag;

  if (n <= 0) {
    ans = 0;
  }
  else if (n == 1) {
    ans = 1;
  }
  else {
    smax = 0.;
    for (i = 1; i <= n; i++) {
      xmag = fabs(SX(i));
      if (smax < xmag) {
        smax = xmag;
        ans  = i;
      }
    }
  }

  return ans;
}

/*============================= end of c_isamax() ========================*/
