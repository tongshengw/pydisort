// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgbfa() ================================*/

/*
    Factors a real band matrix by elimination.
    Revision date:  8/1/82
    Author:  Moler, C. B. (U. of New Mexico)
    c_sgbfa is usually called by c_sgbco, but it can be called
    directly with a saving in time if rcond is not needed.

    Inputs:  same as c_sgbco
    Outputs:
        abd,ipvt    same as c_sgbco
        info    int,
                = 0  normal value.
                = k  if  u(k,k) == 0.  This is not an error
                     condition for this subroutine, but it does
                     indicate that sgbsl will divide by zero if
                     called.  Use  rcond  in c_sgbco for a reliable
                     indication of singularity.
    (see c_sgbco for description of band storage mode)

    NOTE: using memset() to zero columns in abd
 ----------------------------------------------------------------*/

void c_sgbfa(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             int    *info)
{
  register int
    i0,j,j0,j1,ju,jz,k,kp1,l,lm,m,mm,nm1;
  double
    t;

  m     = ml+mu+1;
  *info = 0;
  /*
   * zero initial fill-in columns
   */
  j0 = mu+2;
  j1 = IMIN(n,m)-1;
  for (jz = j0; jz <= j1; jz++) {
    i0 = m+1-jz;
    memset(&ABD(i0,jz),0,(ml-i0+1)*sizeof(double));
  }
  jz = j1;
  ju = 0;

  /*
   * Gaussian elimination with partial pivoting
   */
  nm1 = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
   /*
    * zero next fill-in column
    */
    jz++;
    if (jz <= n) {
      memset(&ABD(1,jz),0,ml*sizeof(double));
    }
    /*
     * find L = pivot index
     */
    lm      = IMIN(ml,n-k);
    l       = c_isamax(lm+1,&ABD(m,k))+m-1;
    IPVT(k) = l+k-m;
    if (ABD(l,k) == 0.) {
     /*
      * zero pivot implies this column already triangularized
      */
      *info = k;
    }
    else {
      /*
       * interchange if necessary
       */
      if (l != m) {
        t        = ABD(l,k);
        ABD(l,k) = ABD(m,k);
        ABD(m,k) = t;
      }
      /*
       * compute multipliers
       */
      t = -1./ABD(m,k);
      c_sscal(lm,t,&ABD(m+1,k));
      /*
       * row elimination with column indexing
       */
      ju = IMIN(IMAX(ju,mu+IPVT(k)),n);
      mm = m;
      for (j = kp1; j <= ju; j++) {
        l--;
        mm--;
        t = ABD(l,j);
        if (l != mm) {
          ABD(l,j)  = ABD(mm,j);
          ABD(mm,j) = t;
        }
        c_saxpy(lm,t,&ABD(m+1,k),&ABD(mm+1,j));
      }
    }
  }
  IPVT(n) = n;
  if (ABD(m,n) == 0.) {
    *info = n;
  }

  return;
}

/*============================= end of c_sgbfa() =========================*/
