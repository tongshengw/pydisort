// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgesl() ================================*/

/*
  Solves the real system
     A*X = B  or  transpose(A)*X = B
  using the factors computed by sgeco or sgefa.
  Revision date:  8/1/82
  Author:  Moler, C. B. (Univ. of New Mexico)

     Inputs:
        a       double(lda, n), the output from sgeco or sgefa.
        lda     int, the leading dimension of the array  A
        n       int, the order of the matrix  A
        ipvt    int(n), the pivot vector from sgeco or sgefa.
        b       double(n), the right hand side vector.
        job     int,
                = 0         to solve  A*X = B ,
                = nonzero   to solve  transpose(A)*X = B

     Outputs:
        b       the solution vector x

     Error condition:
        A division by zero will occur if the input factor contains a
        zero on the diagonal. Technically, this indicates singularity,
        but it is often caused by improper arguments or improper setting
        of lda. It will not occur if the subroutines are called correctly
        and if sgeco has set rcond > 0. or sgefa has set info = 0 .
     To compute  inverse(a)*c where c is a matrix with p columns
           c_sgeco(a,lda,n,ipvt,rcond,z);
           if (rcond is too small) ...
           for (j = 1; j <= p; j++) {
             c_sgesl(a,lda,n,ipvt,c(1,j),0);
           }
 ---------------------------------------------------------------------*/

void c_sgesl(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *b,
             int     job)
{
  register int
    k,kb,l,nm1;
  double
    t;

  nm1 = n-1;
  if (job == 0) {
    /*
     * solve  A*X = B; first solve L*Y = B
     */
    for (k = 1; k <= nm1; k++) {
      l = IPVT(k);
      t = B(l);
      if (l != k) {
        B(l) = B(k);
        B(k) = t;
      }
      c_saxpy(n-k,t,&A(k+1,k),&B(k+1));
    }
    /*
     * now solve  U*X = Y
     */
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      B(k) /= A(k,k);
      t     = -B(k);
      c_saxpy(k-1,t,&A(1,k),&B(1));
    }
  }
  else {
    /*
     * solve trans(A)*X = B; first solve trans(U)*Y = B
     */
    for (k = 1; k <= n; k++) {
      t    = c_sdot(k-1,&A(1,k),&B(1));
      B(k) = (B(k)-t)/A(k,k);
    }
    /*
     * now solve  trans(l)*x = y
     */
    for (kb = 1; kb <= nm1; kb++) {
      k = n-kb;
      B(k) += c_sdot(n-k,&A(k+1,k),&B(k+1));
      l     = IPVT(k);
      if (l != k) {
        t    = B(l);
        B(l) = B(k);
        B(k) = t;
      }
    }
  }

  return;
}

/*============================= end of c_sgesl() =========================*/
