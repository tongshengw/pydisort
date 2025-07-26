// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgeco() ================================*/

/*
   Factors a real matrix by Gaussian elimination
   and estimates the condition of the matrix.
   Revision date:  8/1/82
   Author:  Moler, C. B. (Univ. of New Mexico)
   If rcond is not needed, sgefa is slightly faster.
   To solve  A*X = B, follow sgeco by sgesl.

     Inputs:
        a       double(lda, n), the matrix to be factored.
        lda     int, the leading dimension of the array a.
        n       int, the order of the matrix a.

     Outputs:
        a       an upper triangular matrix and the multipliers
                which were used to obtain it.
                The factorization can be written  A = L*U , where
                L  is a product of permutation and unit lower
                triangular matrices and U is upper triangular.
        ipvt    int(n), an integer vector of pivot indices.
        rcond   double, an estimate of the reciprocal condition of a.
                For the system A*X = B, relative perturbations
                in A and B of size epsilon may cause relative
                perturbations in X of size epsilon/rcond.
                If rcond is so small that the logical expression
                  1.+rcond == 1.
                is true, then A may be singular to working precision.
                In particular, rcond is zero if exact singularity
                is detected or the estimate underflows.
        z       double(n), a work vector whose contents are usually
                unimportant. If A is close to a singular matrix, then z
                is an approximate null vector in the sense that
                norm(A*Z) = rcond*norm(A)*norm(Z) .
 ------------------------------------------------------------------*/

void c_sgeco(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *rcond,
             double *z)
{
  int
    info;
  register int
    j,k,kb,kp1,l;
  double
    anorm,ek,s,sm,t,wk,wkm,ynorm;

  /*
   * compute 1-norm of A
   */
  anorm = 0.;
  for (j = 1; j <= n; j++) {
    anorm = MAX(anorm,c_sasum(n,&A(1,j)));
  }

  /*
   * factor
   */
  c_sgefa(a,lda,n,ipvt,&info);

  /*
   * rcond = 1/(norm(A)*(estimate of norm(inverse(A)))).
   * estimate = norm(Z)/norm(Y) where A*Z = Y and trans(A)*Y = E.
   * trans(A) is the transpose of A. The components of E are
   * chosen to cause maximum local growth in the elements of W where
   * trans(U)*W = E.  The vectors are frequently rescaled to avoid overflow.
   * solve trans(U)*W = E
   */
  ek = 1.;
  memset(z,0,n*sizeof(double));

  for (k = 1; k <= n; k++) {
    if (Z(k) != 0.) {
      ek = F77_SIGN(ek,-Z(k));
    }
    if (fabs(ek-Z(k)) > fabs(A(k,k))) {
      s = fabs(A(k,k))/fabs(ek-Z(k));
      c_sscal(n,s,z);
      ek *= s;
    }
    wk  =  ek-Z(k);
    wkm = -ek-Z(k);
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (A(k,k) != 0.) {
      wk  /= A(k,k);
      wkm /= A(k,k);
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    if (kp1 <= n) {
      for (j = kp1; j <= n; j++) {
        sm   += fabs(Z(j)+wkm*A(k,j));
        Z(j) += wk*A(k,j);
        s    += fabs(Z(j));
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        for (j = kp1; j <= n; j++) {
          Z(j) += t*A(k,j);
        }
      }
    }
    Z(k) = wk;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve trans(L)*Y = W
   */
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (k < n) {
      Z(k) += c_sdot(n-k,&A(k+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
    }
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
  }
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve L*V = Y
   */
  ynorm = 1.;
  for (k = 1; k <= n; k++) {
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
    if (k < n) {
      c_saxpy(n-k,t,&A(k+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
  }
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve U*Z = V
   */
  ynorm *= s;
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(Z(k)) > fabs(A(k,k))) {
      s = fabs(A(k,k))/fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
    if (A(k,k) != 0.) {
      Z(k) /= A(k,k);
    }
    else {
      Z(k) = 1.;
    }
    t = -Z(k);
    c_saxpy(k-1,t,&A(1,k),&Z(1));
  }
  /*
   * make znorm = 1.0
   */
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  ynorm *= s;
  if (anorm != 0.) {
    *rcond = ynorm/anorm;
  }
  else {
    *rcond = 0.;
  }

  return;
}

/*============================= end of c_sgeco() =========================*/
