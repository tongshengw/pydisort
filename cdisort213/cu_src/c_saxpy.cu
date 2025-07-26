// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_saxpy() ================================*/

/*
  y = a*x + y  (x, y = vectors, a = scalar)

  INPUT--
        n   Number of elements in input vectors x and y
       sa   Scalar multiplier a
       sx   Array containing vector x
       sy   Array containing vector Y

 OUTPUT--
       sy   For i = 1 to n, overwrite  SY(i) with sa*SX(i)+SY(i)

  NOTE: Fortran inputs incx, incy removed because they are not used
        by disort or twostr
 ------------------------------------------------------------*/

void c_saxpy(int     n,
             double  sa,
             double *sx,
             double *sy)
{
  register int
    i,m;

  if (n <= 0 || sa == 0.) {
    return;
  }

  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      SY(i) += sa*SX(i);
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    SY(i  ) += sa*SX(i  );
    SY(i+1) += sa*SX(i+1);
    SY(i+2) += sa*SX(i+2);
    SY(i+3) += sa*SX(i+3);
  }

  return;
}

/*============================= end of c_saxpy() =========================*/
