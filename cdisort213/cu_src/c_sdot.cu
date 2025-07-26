// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sdot() =================================*/

/*
  Dot product of vectors x and y

  INPUT--
        n  Number of elements in input vectors x and y
       sx  Array containing vector x
       sy  Array containing vector y

 OUTPUT--
      ans  Sum for i = 1 to n of  SX(i)*SY(i),

  NOTE: Fortran input arguments incx, incy removed because they
        are not used in disort or twostr
 ------------------------------------------------------------------*/

double c_sdot(int     n,
              double *sx,
              double *sy)
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
      ans += SX(i)*SY(i);
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    ans += SX(i  )*SY(i  )
          +SX(i+1)*SY(i+1)
          +SX(i+2)*SY(i+2)
          +SX(i+3)*SY(i+3);
  }

  return ans;
}

/*============================= end of c_sdot() ==========================*/
