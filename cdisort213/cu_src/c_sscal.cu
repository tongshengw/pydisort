// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sscal() ================================*/

/*
  Multiply vector sx by scalar sa

  INPUT--  n  Number of elements in vector
          sa  Scale factor
          sx  Array, length n, containing vector

 OUTPUT-- sx  Replace SX(i) with sa*SX(i) for i = i to n

 NOTE: Fortran input argument incx removed since it is not used
       in disort or twostr

 ---------------------------------------------------------------------*/

void c_sscal(int    n,
             double  sa,
             double *sx)
{
  register int
    i,m;

  if (n <= 0) {
    return;
  }
  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      SX(i) *= sa;
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    SX(i  ) *= sa;
    SX(i+1) *= sa;
    SX(i+2) *= sa;
    SX(i+3) *= sa;
  }

  return;
}

/*============================= end of c_sscal() =========================*/
