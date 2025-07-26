// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sasum() ================================*/

/*
  Input--   n     Number of elements in vector to be summed
            sx    array, length n, containing vector

  OUTPUT--  ans   Sum from i = 1 to n of fabs(SX(i))

  NOTE: Fortran input incx removed because it is not used by
        disort or twostr
 ----------------------------------------------------------*/

double c_sasum(int     n,
             double *sx)
{
  register int
    i,m;
  double
    ans;

  ans = 0.;
  if (n <= 0) {
    return ans;
  }

  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      ans += fabs(SX(i));
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    ans += fabs(SX(i  ))
          +fabs(SX(i+1))
          +fabs(SX(i+2))
          +fabs(SX(i+3));
  }

  return ans;
}

/*============================= end of c_sasum() =========================*/
