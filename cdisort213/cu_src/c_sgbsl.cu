// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgbsl() ================================*/

/*
    Solves the real band system
       A * X = B  or  transpose(A) * X = B
    using the factors computed by sgbco or sgbfa.
    Revision date:  8/1/82
    Author:  Moler, C. B. (Univ. of New Mexico)

    Inputs:
        abd     double(lda, n), the output from sgbco or sgbfa.
        lda     int, the leading dimension of the array abd.
        n       int, the order of the original matrix.
        ml      int, number of diagonals below the main diagonal.
        mu      int, number of diagonals above the main diagonal.
        ipvt    int(n), the pivot vector from sgbco or sgbfa.
        b       double(n), the right hand side vector.
        job     int,
                = 0         to solve  A*X = B ,
                = nonzero   to solve  transpose(A)*X = B

     Outputs:
        b       the solution vector  X

     Error condition:
        A division by zero will occur if the input factor contains a
        zero on the diagonal.  Technically, this indicates singularity,
        but it is often caused by improper arguments or improper
        setting of lda.  It will not occur if the subroutines are
        called correctly and if c_sgbco has set rcond > 0.0
        or sgbfa has set info = 0 .
     To compute  inverse(a)*c  where c is a matrix
     with p columns
      c_sgbco(abd,lda,n,ml,mu,ipvt,&rcond,z)
      if (rcond is too small) ...
        for (j = 1; j <= p; j++) {
          c_sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
        }
 --------------------------------------------------------*/

void c_sgbsl(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *b,
             int     job)
{
  register int
    k,kb,l,la,lb,lm,m,nm1;
  double
    t;

  m   = mu+ml+1;
  nm1 = n-1;
  if (job == 0) {
   /*
    * solve  A*X = B;  first solve L*Y = B
    */
    if (ml != 0) {
      for (k = 1; k <= nm1; k++) {
        lm = IMIN(ml,n-k);
        l  = IPVT(k);
        t  = B(l);
        if (l != k) {
          B(l) = B(k);
          B(k) = t;
        }
        c_saxpy(lm,t,&ABD(m+1,k),&B(k+1));
      }
    }
    /*
     * now solve  U*X = Y
     */
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      B(k) /= ABD(m,k);
      lm    = IMIN(k,m)-1;
      la    = m-lm;
      lb    = k-lm;
      t     = -B(k);
      c_saxpy(lm,t,&ABD(la,k),&B(lb));
    }
  }
  else {
    /*
     * solve  trans(A)*X = B; first solve trans(U)*Y = B
     */
    for (k = 1; k <= n; k++) {
      lm   = IMIN(k,m)-1;
      la   = m-lm;
      lb   = k-lm;
      t    = c_sdot(lm,&ABD(la,k),&B(lb));
      B(k) = (B(k)-t)/ABD(m,k);
    }
    /*
     * now solve trans(L)*X = Y
     */
    if (ml != 0) {
      for (kb = 1; kb <= nm1; kb++) {
        k     = n-kb;
        lm    = IMIN(ml,n-k);
        B(k) += c_sdot(lm,&ABD(m+1,k),&B(k+1));
        l     = IPVT(k);
        if (l != k) {
          t    = B(l);
          B(l) = B(k);
          B(k) = t;
        }
      }
    }
  }

  return;
}

/*============================= end of c_sgbsl() =========================*/
