// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgefa() ================================*/

/*
   Factors a real matrix by Gaussian elimination.
   Revision date:  8/1/82
   Author:  Moler, C. B. (Univ. of New Mexico)
   c_sgefa is usually called by c_sgeco, but it can be called directly with a
   saving in time if rcond is not needed.
   (time for c_sgeco) = (1+9/n)*(time for c_sgefa).

   Inputs:  same as c_sgeco

   Outputs:
        a,ipvt  same as c_sgeco
        info    int,
                = 0  normal value.
                = k  if  u(k,k) = 0.  This is not an error condition for
                     this subroutine, but it does indicate that c_sgesl or
                     c_sgedi will divide by zero if called.  Use rcond in
                     c_sgeco for a reliable indication of singularity.
 ---------------------------------------------------------------------*/

void c_sgefa(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             int    *info)
{
  register int
    j,k,kp1,l,nm1;
  double
    t;

  /*
   * Gaussian elimination with partial pivoting
   */
  *info = 0;
  nm1   = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
    /*
     * find L = pivot index
     */
    l       = c_isamax(n-k+1,&A(k,k))+k-1;
    IPVT(k) = l;
    if (A(l,k) == 0.) {
      /*
       * zero pivot implies this column already triangularized
       */
      *info = k;
    }
    else {
      /*
       * interchange if necessary
       */
      if (l != k) {
        t      = A(l,k);
        A(l,k) = A(k,k);
        A(k,k) = t;
      }
      /*
       * compute multipliers
       */
      t = -1./A(k,k);
      c_sscal(n-k,t,&A(k+1,k));
      /*
       * row elimination with column indexing
       */
      for (j = kp1; j <= n; j++) {
        t = A(l,j);
        if (l != k) {
          A(l,j) = A(k,j);
          A(k,j) = t;
        }
        c_saxpy(n-k,t,&A(k+1,k),&A(k+1,j));
      }
    }
  }
  IPVT(n) = n;
  if (A(n,n) == 0.) {
    *info = n;
  }

  return;
}

/*============================= end of c_sgefa() =========================*/
